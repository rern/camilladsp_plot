import cmath
import math
import math
from .cooley_tukey import fft
from .audiofileread import read_coeffs


def unwrap_phase(values, threshold=150.0):
    offset = 0
    prevdiff = 0.0
    unwrapped = [0.0] * len(values)
    if len(values) > 0:
        unwrapped[0] = values[0]
        for n in range(1, len(values)):
            guess = values[n - 1] + prevdiff
            diff = values[n] - guess
            if diff > threshold:
                offset -= 1
                jumped = True
            elif diff < -threshold:
                offset += 1
                jumped = True
            else:
                jumped = False
            unwrapped[n] = values[n] + 2 * 180.0 * offset
            if not jumped:
                prevdiff = unwrapped[n] - unwrapped[n - 1]
    return unwrapped


def calc_groupdelay(freq, phase):
    if len(freq) < 2:
        return [], []
    phase = unwrap_phase(phase)
    freq_new = []
    groupdelay = []
    for n in range(1, len(freq)):
        dw = (freq[n] - freq[n - 1]) * 2 * math.pi
        f = (freq[n - 1] + freq[n]) / 2.0
        freq_new.append(f)
        dp = (phase[n] - phase[n - 1]) / 180.0 * math.pi
        delay = -1000.0 * dp / dw
        groupdelay.append(delay)
    return freq_new, groupdelay


class BaseFilter(object):
    def __init__(self):
        pass

    def complex_gain(self, f, remove_delay=False):
        A = [1.0] * len(f)
        return f, A

    def gain_and_phase(self, f, remove_delay=False):
        _f, Avec = self.complex_gain(f, remove_delay=remove_delay)
        gain = [20 * math.log10(abs(A) + 1.0e-15) for A in Avec]
        phase = [180 / math.pi * cmath.phase(A) for A in Avec]
        return f, gain, phase

    def is_stable(self):
        return True


class Conv(object):

    def __init__(self, conf, fs):
        if not conf:
            conf = {"values": [1.0]}
        if "filename" in conf:
            values = read_coeffs(conf)
        else:
            values = conf["values"]
        self.impulse = values
        self.fs = fs

    def find_peak(self):
        (_, idx) = max((abs(val), idx) for idx, val in enumerate(self.impulse))
        return idx

    def complex_gain(self, f, remove_delay=False):
        impulselen = len(self.impulse)
        npoints = 2 ** (math.ceil(math.log2(impulselen)))
        if npoints < 1024:
            npoints = 1024
        impulse = list(self.impulse)
        padding = [0.0 for _n in range(npoints - impulselen)]
        impulse.extend(padding)
        impfft = fft(impulse)
        f_fft = [self.fs * n / (npoints) for n in range(int(npoints / 2))]
        cut = impfft[0 : int(npoints / 2)]
        if remove_delay:
            maxidx = self.find_peak()
            cut = [
                val * cmath.exp(1j * 1.0 / (npoints / 2) * math.pi * idx * maxidx)
                for idx, val in enumerate(cut)
            ]
        if f is not None:
            interpolated = self.interpolate_polar(cut, f_fft, f)
            return f, interpolated
        return f_fft, cut

    def interpolate(self, y, xold, xnew):
        idx = 0
        ynew = []

        for x in xnew:
            idx = len(y) * x / xold[-1]
            i1 = int(math.floor(idx))
            i2 = i1 + 1
            if i1 >= (len(y)):
                i1 = len(y) - 1
            if i2 >= (len(y)):
                i2 = i1
            fract = idx - i1
            newval = (1 - fract) * y[i1] + fract * y[i2]
            ynew.append(newval)
        return ynew

    def interpolate_polar(self, y, xold, xnew):
        y_magn = [abs(yval) for yval in y]
        y_ang = [180.0 / math.pi * cmath.phase(yval) for yval in y]
        y_ang = [
            math.pi * yval / 180.0 for yval in unwrap_phase(y_ang, threshold=270.0)
        ]
        y_magn_interp = self.interpolate(y_magn, xold, xnew)
        y_ang_interp = self.interpolate(y_ang, xold, xnew)
        return [cmath.rect(r, phi) for (r, phi) in zip(y_magn_interp, y_ang_interp)]

    def gain_and_phase(self, f, remove_delay=False):
        f_fft, Avec = self.complex_gain(None, remove_delay=remove_delay)
        interpolated = self.interpolate_polar(Avec, f_fft, f)
        gain = [20.0 * math.log10(abs(A) + 1.0e-15) for A in interpolated]
        phase = [180.0 / math.pi * cmath.phase(A) for A in interpolated]
        return f, gain, phase

    def get_impulse(self):
        t = [n / self.fs for n in range(len(self.impulse))]
        return t, self.impulse


class DiffEq(BaseFilter):
    def __init__(self, conf, fs):
        self.fs = fs
        self.a = conf["a"]
        self.b = conf["b"]
        if len(self.a) == 0:
            self.a = [1.0]
        if len(self.b) == 0:
            self.b = [1.0]

    def complex_gain(self, freq, remove_delay=False):
        zvec = [cmath.exp(1j * 2 * math.pi * f / self.fs) for f in freq]
        A1 = [0.0 for n in range(len(freq))]
        for n, bn in enumerate(self.b):
            A1 = [a1 + bn * z ** (-n) for a1, z in zip(A1, zvec)]
        A2 = [0.0 for n in range(len(freq))]
        for n, an in enumerate(self.a):
            A2 = [a2 + an * z ** (-n) for a2, z in zip(A2, zvec)]
        A = [a1 / a2 for (a1, a2) in zip(A1, A2)]
        return freq, A

    def is_stable(self):
        # TODO
        return None


class Delay(BaseFilter):
    def __init__(self, conf, fs):
        self.fs = fs
        unit = conf.get("unit", "ms")
        if unit is None:
            unit = "ms"
        if unit == "ms":
            self.delay_samples = conf["delay"] / 1000.0 * fs
        elif unit == "mm":
            self.delay_samples = conf["delay"] / 1000.0 * fs / 343.0
        elif unit == "samples":
            self.delay_samples = conf["delay"]
        else:
            raise RuntimeError(f"Unknown unit {unit}")

        self.subsample = conf.get("subsample", False) == True
        if self.subsample:
            self.delay_full_samples = math.floor(self.delay_samples)
            self.fraction = self.delay_samples - self.delay_full_samples
            self.a1 = 1.0 - self.fraction
            self.a2 = 0.0
            self.b0 = 1.0 - self.fraction
            self.b1 = 1.0
            self.b2 = 0.0
        else:
            self.delay_full_samples = round(self.delay_samples)

    def complex_gain(self, freq, remove_delay=False):
        zvec = [cmath.exp(1j * 2 * math.pi * f / self.fs) for f in freq]
        if self.subsample:
            A = [
                (
                    (self.b0 + self.b1 * z ** (-1) + self.b2 * z ** (-2))
                    / (1.0 + self.a1 * z ** (-1) + self.a2 * z ** (-2))
                )
                for z in zvec
            ]
        else:
            A = [1.0 for _z in zvec]
        if not remove_delay:
            delay_s = self.delay_full_samples / self.fs
            A = [
                val * cmath.exp(-1j * 2.0 * math.pi * f * delay_s)
                for val, f in zip(A, freq)
            ]
        return freq, A

    def is_stable(self):
        # TODO
        return None


class Gain(BaseFilter):
    def __init__(self, conf):
        self.gain = conf["gain"]
        self.inverted = conf["inverted"] == True
        self.scale = conf["scale"]
        if self.scale is None:
            self.scale = "dB"

    def complex_gain(self, f, remove_delay=False):
        sign = -1.0 if self.inverted else 1.0
        if self.scale == "dB":
            gain = 10.0 ** (self.gain / 20.0) * sign
        else:
            gain = self.gain * sign
        A = [gain for n in range(len(f))]
        return f, A


class BiquadCombo(BaseFilter):
    def Butterw_q(self, order):
        odd = order % 2 > 0
        n_so = math.floor(order / 2.0)
        qvalues = []
        for n in range(0, n_so):
            q = 1 / (2.0 * math.sin((math.pi / order) * (n + 1 / 2)))
            qvalues.append(q)
        if odd:
            qvalues.append(-1.0)
        return qvalues

    def __init__(self, conf, fs):
        self.ftype = conf["type"]
        if self.ftype in [
            "LinkwitzRileyHighpass",
            "LinkwitzRileyLowpass",
            "ButterworthHighpass",
            "ButterworthHighpass",
            "ButterworthLowpass",
        ]:
            self.order = conf["order"]
            self.freq = conf["freq"]
            self.fs = fs
            if self.ftype == "LinkwitzRileyHighpass":
                # qvalues = self.LRtable[self.order]
                q_temp = self.Butterw_q(self.order / 2)
                if (self.order / 2) % 2 > 0:
                    q_temp = q_temp[0:-1]
                    qvalues = q_temp + q_temp + [0.5]
                else:
                    qvalues = q_temp + q_temp
                type_so = "Highpass"
                type_fo = "HighpassFO"

            elif self.ftype == "LinkwitzRileyLowpass":
                q_temp = self.Butterw_q(self.order / 2)
                if (self.order / 2) % 2 > 0:
                    q_temp = q_temp[0:-1]
                    qvalues = q_temp + q_temp + [0.5]
                else:
                    qvalues = q_temp + q_temp
                type_so = "Lowpass"
                type_fo = "LowpassFO"
            elif self.ftype == "ButterworthHighpass":
                qvalues = self.Butterw_q(self.order)
                type_so = "Highpass"
                type_fo = "HighpassFO"
            elif self.ftype == "ButterworthLowpass":
                qvalues = self.Butterw_q(self.order)
                type_so = "Lowpass"
                type_fo = "LowpassFO"
            self.biquads = []
            for q in qvalues:
                if q >= 0:
                    bqconf = {"freq": self.freq, "q": q, "type": type_so}
                else:
                    bqconf = {"freq": self.freq, "type": type_fo}
                self.biquads.append(Biquad(bqconf, self.fs))
        elif self.ftype == "FivePointPeq":
            lsconf = Biquad(
                {
                    "freq": conf["fls"],
                    "q": conf["qls"],
                    "gain": conf["gls"],
                    "type": "Lowshelf",
                },
                fs,
            )
            hsconf = Biquad(
                {
                    "freq": conf["fhs"],
                    "q": conf["qhs"],
                    "gain": conf["ghs"],
                    "type": "Highshelf",
                },
                fs,
            )
            p1conf = Biquad(
                {
                    "freq": conf["fp1"],
                    "q": conf["qp1"],
                    "gain": conf["gp1"],
                    "type": "Peaking",
                },
                fs,
            )
            p2conf = Biquad(
                {
                    "freq": conf["fp2"],
                    "q": conf["qp2"],
                    "gain": conf["gp2"],
                    "type": "Peaking",
                },
                fs,
            )
            p3conf = Biquad(
                {
                    "freq": conf["fp3"],
                    "q": conf["qp3"],
                    "gain": conf["gp3"],
                    "type": "Peaking",
                },
                fs,
            )
            self.biquads = [lsconf, p1conf, p2conf, p3conf, hsconf]
        elif self.ftype == "GraphicEqualizer":
            bands = len(conf["gains"])
            f_min = conf["freq_min"] if conf["freq_min"] else 20.0
            f_max = conf["freq_max"] if conf["freq_max"] else 20000.0
            f_min_log = math.log2(f_min)
            f_max_log = math.log2(f_max)
            self.biquads = []
            bw = (f_max_log - f_min_log) / bands
            for band, gain in enumerate(conf["gains"]):
                if math.fabs(gain) > 0.01:
                    freq_log = f_min_log + (band + 0.5) * bw
                    freq = 2.0**freq_log
                    filt = Biquad(
                        {
                            "freq": freq,
                            "bandwidth": bw,
                            "gain": gain,
                            "type": "Peaking",
                        },
                        fs,
                    )
                    self.biquads.append(filt)
        elif self.ftype == "Tilt":
            gain_low = -conf["gain"] / 2.0
            gain_high = conf["gain"] / 2.0
            lsconf = Biquad(
                {"freq": 110.0, "q": 0.35, "gain": gain_low, "type": "Lowshelf"}, fs
            )
            hsconf = Biquad(
                {"freq": 3500.0, "q": 0.35, "gain": gain_high, "type": "Highshelf"}, fs
            )
            self.biquads = [lsconf, hsconf]

    def is_stable(self):
        # TODO
        return None

    def complex_gain(self, freq, remove_delay=False):
        A = [1.0 for n in range(len(freq))]
        for bq in self.biquads:
            _f, Atemp = bq.complex_gain(freq)
            A = [a * atemp for (a, atemp) in zip(A, Atemp)]
        return freq, A


class Loudness(BaseFilter):
    def __init__(self, conf, fs, volume):
        rel_vol = volume - conf["reference_level"]
        conf["low_boost"]
        conf["attenuate_mid"]
        rel_boost = -rel_vol / 20.0
        if rel_boost > 1.0:
            rel_boost = 1.0
        elif rel_boost < 0.0:
            rel_boost = 0.0
        high_boost = rel_boost * conf["high_boost"]
        low_boost = rel_boost * conf["low_boost"]
        if conf["attenuate_mid"]:
            max_gain = max(high_boost, low_boost)
            self.mid_gain = 10.0 ** (-max_gain / 20.0)
        else:
            self.mid_gain = 1.0

        lsconf = Biquad(
            {"freq": 70.0, "slope": 12.0, "gain": low_boost, "type": "Lowshelf"}, fs
        )
        hsconf = Biquad(
            {"freq": 3500.0, "slope": 12.0, "gain": high_boost, "type": "Highshelf"}, fs
        )
        self.biquads = [lsconf, hsconf]

    def complex_gain(self, freq, remove_delay=False):
        A = [self.mid_gain for n in range(len(freq))]
        for bq in self.biquads:
            _f, Atemp = bq.complex_gain(freq)
            A = [a * atemp for (a, atemp) in zip(A, Atemp)]
        return freq, A


class Biquad(BaseFilter):
    def __init__(self, conf, fs):
        ftype = conf["type"]
        if ftype == "Free":
            a0 = 1.0
            a1 = conf["a1"]
            a2 = conf["a2"]
            b0 = conf["b0"]
            b1 = conf["b1"]
            b2 = conf["b2"]
        if ftype == "Highpass":
            freq = conf["freq"]
            q = conf["q"]
            omega = 2.0 * math.pi * freq / fs
            sn = math.sin(omega)
            cs = math.cos(omega)
            alpha = sn / (2.0 * q)
            b0 = (1.0 + cs) / 2.0
            b1 = -(1.0 + cs)
            b2 = (1.0 + cs) / 2.0
            a0 = 1.0 + alpha
            a1 = -2.0 * cs
            a2 = 1.0 - alpha
        elif ftype == "Lowpass":
            freq = conf["freq"]
            q = conf["q"]
            omega = 2.0 * math.pi * freq / fs
            sn = math.sin(omega)
            cs = math.cos(omega)
            alpha = sn / (2.0 * q)
            b0 = (1.0 - cs) / 2.0
            b1 = 1.0 - cs
            b2 = (1.0 - cs) / 2.0
            a0 = 1.0 + alpha
            a1 = -2.0 * cs
            a2 = 1.0 - alpha
        elif ftype == "Peaking":
            freq = conf["freq"]
            gain = conf["gain"]
            omega = 2.0 * math.pi * freq / fs
            sn = math.sin(omega)
            cs = math.cos(omega)
            ampl = 10.0 ** (gain / 40.0)
            if "q" in conf:
                q = conf["q"]
                alpha = sn / (2.0 * q)
            else:
                bandwidth = conf["bandwidth"]
                alpha = sn * math.sinh(math.log(2.0) / 2.0 * bandwidth * omega / sn)
            b0 = 1.0 + (alpha * ampl)
            b1 = -2.0 * cs
            b2 = 1.0 - (alpha * ampl)
            a0 = 1.0 + (alpha / ampl)
            a1 = -2.0 * cs
            a2 = 1.0 - (alpha / ampl)
        elif ftype == "HighshelfFO":
            freq = conf["freq"]
            gain = conf["gain"]
            omega = 2.0 * math.pi * freq / fs
            ampl = 10.0 ** (gain / 40.0)
            tn = math.tan(omega / 2)
            b0 = ampl * tn + ampl**2
            b1 = ampl * tn - ampl**2
            b2 = 0.0
            a0 = ampl * tn + 1
            a1 = ampl * tn - 1
            a2 = 0.0
        elif ftype == "Highshelf":
            freq = conf["freq"]
            gain = conf["gain"]
            omega = 2.0 * math.pi * freq / fs
            ampl = 10.0 ** (gain / 40.0)
            sn = math.sin(omega)
            cs = math.cos(omega)
            if "slope" in conf:
                slope = conf["slope"]
                alpha = (
                    sn
                    / 2.0
                    * math.sqrt(
                        (ampl + 1.0 / ampl) * (1.0 / (slope / 12.0) - 1.0) + 2.0
                    )
                )
                beta = 2.0 * math.sqrt(ampl) * alpha
            else:
                q = conf["q"]
                beta = sn * math.sqrt(ampl) / q
            b0 = ampl * ((ampl + 1.0) + (ampl - 1.0) * cs + beta)
            b1 = -2.0 * ampl * ((ampl - 1.0) + (ampl + 1.0) * cs)
            b2 = ampl * ((ampl + 1.0) + (ampl - 1.0) * cs - beta)
            a0 = (ampl + 1.0) - (ampl - 1.0) * cs + beta
            a1 = 2.0 * ((ampl - 1.0) - (ampl + 1.0) * cs)
            a2 = (ampl + 1.0) - (ampl - 1.0) * cs - beta
        elif ftype == "LowshelfFO":
            freq = conf["freq"]
            gain = conf["gain"]
            omega = 2.0 * math.pi * freq / fs
            ampl = 10.0 ** (gain / 40.0)
            tn = math.tan(omega / 2)
            b0 = ampl**2 * tn + ampl
            b1 = ampl**2 * tn - ampl
            b2 = 0.0
            a0 = tn + ampl
            a1 = tn - ampl
            a2 = 0.0
        elif ftype == "Lowshelf":
            freq = conf["freq"]
            gain = conf["gain"]
            omega = 2.0 * math.pi * freq / fs
            ampl = 10.0 ** (gain / 40.0)
            sn = math.sin(omega)
            cs = math.cos(omega)
            if "slope" in conf:
                slope = conf["slope"]
                alpha = (
                    sn
                    / 2.0
                    * math.sqrt(
                        (ampl + 1.0 / ampl) * (1.0 / (slope / 12.0) - 1.0) + 2.0
                    )
                )
                beta = 2.0 * math.sqrt(ampl) * alpha
            else:
                q = conf["q"]
                beta = sn * math.sqrt(ampl) / q

            b0 = ampl * ((ampl + 1.0) - (ampl - 1.0) * cs + beta)
            b1 = 2.0 * ampl * ((ampl - 1.0) - (ampl + 1.0) * cs)
            b2 = ampl * ((ampl + 1.0) - (ampl - 1.0) * cs - beta)
            a0 = (ampl + 1.0) + (ampl - 1.0) * cs + beta
            a1 = -2.0 * ((ampl - 1.0) + (ampl + 1.0) * cs)
            a2 = (ampl + 1.0) + (ampl - 1.0) * cs - beta
        elif ftype == "LowpassFO":
            freq = conf["freq"]
            omega = 2.0 * math.pi * freq / fs
            k = math.tan(omega / 2.0)
            alpha = 1 + k
            a0 = 1.0
            a1 = -((1 - k) / alpha)
            a2 = 0.0
            b0 = k / alpha
            b1 = k / alpha
            b2 = 0
        elif ftype == "HighpassFO":
            freq = conf["freq"]
            omega = 2.0 * math.pi * freq / fs
            k = math.tan(omega / 2.0)
            alpha = 1 + k
            a0 = 1.0
            a1 = -((1 - k) / alpha)
            a2 = 0.0
            b0 = 1.0 / alpha
            b1 = -1.0 / alpha
            b2 = 0
        elif ftype == "Notch":
            freq = conf["freq"]
            omega = 2.0 * math.pi * freq / fs
            sn = math.sin(omega)
            cs = math.cos(omega)
            if "q" in conf:
                q = conf["q"]
                alpha = sn / (2.0 * q)
            else:
                bandwidth = conf["bandwidth"]
                alpha = sn * math.sinh(math.log(2.0) / 2.0 * bandwidth * omega / sn)
            b0 = 1.0
            b1 = -2.0 * cs
            b2 = 1.0
            a0 = 1.0 + alpha
            a1 = -2.0 * cs
            a2 = 1.0 - alpha
        elif ftype == "GeneralNotch":
            f_p = conf["freq_pole"]
            f_z = conf["freq_zero"]
            q_p = conf["q_pole"]
            normalize_at_dc = conf["normalize_at_dc"] == True

            # apply pre-warping
            tn_z = math.tan(math.pi * f_z / fs)
            tn_p = math.tan(math.pi * f_p / fs)
            alpha = tn_p / q_p
            tn2_p = tn_p**2
            tn2_z = tn_z**2

            # calculate gain
            if normalize_at_dc:
                gain = tn2_p / tn2_z
            else:
                gain = 1.0

            b0 = gain * (1.0 + tn2_z)
            b1 = -2.0 * gain * (1.0 - tn2_z)
            b2 = gain * (1.0 + tn2_z)
            a0 = 1.0 + alpha + tn2_p
            a1 = -2.0 + 2.0 * tn2_p
            a2 = 1.0 - alpha + tn2_p

        elif ftype == "Bandpass":
            freq = conf["freq"]
            omega = 2.0 * math.pi * freq / fs
            sn = math.sin(omega)
            cs = math.cos(omega)
            if "q" in conf:
                q = conf["q"]
                alpha = sn / (2.0 * q)
            else:
                bandwidth = conf["bandwidth"]
                alpha = sn * math.sinh(math.log(2.0) / 2.0 * bandwidth * omega / sn)
            b0 = alpha
            b1 = 0.0
            b2 = -alpha
            a0 = 1.0 + alpha
            a1 = -2.0 * cs
            a2 = 1.0 - alpha
        elif ftype == "Allpass":
            freq = conf["freq"]
            omega = 2.0 * math.pi * freq / fs
            sn = math.sin(omega)
            cs = math.cos(omega)
            if "q" in conf:
                q = conf["q"]
                alpha = sn / (2.0 * q)
            else:
                bandwidth = conf["bandwidth"]
                alpha = sn * math.sinh(math.log(2.0) / 2.0 * bandwidth * omega / sn)
            b0 = 1.0 - alpha
            b1 = -2.0 * cs
            b2 = 1.0 + alpha
            a0 = 1.0 + alpha
            a1 = -2.0 * cs
            a2 = 1.0 - alpha
        elif ftype == "AllpassFO":
            freq = conf["freq"]
            omega = 2.0 * math.pi * freq / fs
            tn = math.tan(omega / 2.0)
            alpha = (tn + 1.0) / (tn - 1.0)
            b0 = 1.0
            b1 = alpha
            b2 = 0.0
            a0 = alpha
            a1 = 1.0
            a2 = 0.0
        elif ftype == "LinkwitzTransform":
            f0 = conf["freq_act"]
            q0 = conf["q_act"]
            qt = conf["q_target"]
            ft = conf["freq_target"]

            d0i = (2.0 * math.pi * f0) ** 2
            d1i = (2.0 * math.pi * f0) / q0
            c0i = (2.0 * math.pi * ft) ** 2
            c1i = (2.0 * math.pi * ft) / qt
            fc = (ft + f0) / 2.0

            gn = 2 * math.pi * fc / math.tan(math.pi * fc / fs)
            cci = c0i + gn * c1i + gn**2

            b0 = (d0i + gn * d1i + gn**2) / cci
            b1 = 2 * (d0i - gn**2) / cci
            b2 = (d0i - gn * d1i + gn**2) / cci
            a0 = 1.0
            a1 = 2.0 * (c0i - gn**2) / cci
            a2 = (c0i - gn * c1i + gn**2) / cci

        self.fs = fs
        self.a1 = a1 / a0
        self.a2 = a2 / a0
        self.b0 = b0 / a0
        self.b1 = b1 / a0
        self.b2 = b2 / a0

    def complex_gain(self, freq, remove_delay=False):
        zvec = [cmath.exp(1j * 2 * math.pi * f / self.fs) for f in freq]
        A = [
            (
                (self.b0 + self.b1 * z ** (-1) + self.b2 * z ** (-2))
                / (1.0 + self.a1 * z ** (-1) + self.a2 * z ** (-2))
            )
            for z in zvec
        ]
        return freq, A

    def is_stable(self):
        return abs(self.a2) < 1.0 and abs(self.a1) < (self.a2 + 1.0)
