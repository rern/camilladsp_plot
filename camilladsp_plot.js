/*
convert from https://github.com/HEnquist/pycamilladsp-plot/blob/master/camilladsp_plot/
this file - camilladsp_plot.js << filters.py
	- require complex.js - math operations >> serial of complex function
		- Complex( 0, 1 )
			.mul( 2 * Math.PI * f / this.fs )
			.exp()
	- add fft() for Conv complexGain()
camilla.js:
	- graph.filter(), graph.pipeline() << eval_filterconfig.py
camilla_audiofilread.py << audiofilread.py
	- modify for graph.filter() - Conv
*/
function calcGroupDelay( freq, phase ) {
	if ( freq.length < 2 ) {
		return [ [], [] ];
	}
	phase = unwrapPhase( phase );
	var freqNew    = [];
	var groupDelay = [];
	for ( var n = 1; n < freq.length; n++ ) {
		var dw    = ( freq[ n ] - freq[ n - 1 ] ) * 2 * Math.PI;
		var f     = ( freq[ n - 1] + freq[ n ] ) / 2;
		freqNew.push( f );
		var dp    = ( phase[ n ] - phase[ n - 1 ] ) / 180 * Math.PI;
		var delay = -1000 * dp / dw;
		groupDelay.push( delay );
	}
	return [ freqNew, groupDelay ];
}
function mapGain( cgain ) {
	return cgain.map( A => isFinite( A.re ) ? 20 * Math.log10( A.abs() + 1e-15 ) : 0 )
}
function mapPhase( cphase ) {
	return cphase.map( A => Math.atan2( A.im, A.re ) * 180 / Math.PI )
}
function unwrapPhase( values, threshold = 150 ) {
	var offset    = 0;
	var prevdiff  = 0;
	var unwrapped = new Array( values.length ).fill( 0 );
	if ( values.length > 0 ) {
		unwrapped[ 0 ] = values[ 0 ];
		for ( var n = 1; n < values.length; n++ ) {
			var guess = values[ n - 1 ] + prevdiff;
			var diff  = values[ n ] - guess;
			var jumped = false;
			if ( diff > threshold ) {
				offset -= 1;
				jumped  = true;
			} else if ( diff < -threshold ) {
				offset += 1;
				jumped  = true;
			}
			unwrapped[ n ] = values[ n ] + 2 * 180 * offset;
			if ( ! jumped ) {
				prevdiff = unwrapped[ n ] - unwrapped[ n - 1 ];
			}
		}
	}
	return unwrapped;
}

class BaseFilter {
	constructor() {}

	complexGain( freq, removeDelay = false ) {
		var A = new Array( freq.length ).fill( 1 );
		return [ freq, A ];
	}

	gainAndPhase( freq, removeDelay = false, callback ) {
		var [_, Avec]  = this.complexGain(freq, removeDelay);
		var ga         = mapGain( Avec );
		var ph         = mapPhase( Avec );
		var [ fg, gr ] = calcGroupDelay( freq, ph );
		callback( fg, ga, gr, ph );
	}
}
class Biquad extends BaseFilter {
	constructor( conf, fs ) {
		super();
		var ftype = conf.type;
		var a0, a1, a2, b0, b1, b2;
		switch ( ftype ) {
			case 'Allpass':
				var freq  = conf.freq;
				var omega = 2 * Math.PI * freq / fs;
				var sn    = Math.sin( omega );
				var cs    = Math.cos( omega );
				var alpha = 'q' in conf ? sn / ( 2 * conf.q ) : sn * Math.sinh( Math.log( 2 ) / 2 * conf.bandwidth * omega / sn );
				b0        = 1 - alpha;
				b1        = -2 * cs;
				b2        = 1 + alpha;
				a0        = 1 + alpha;
				a1        = -2 * cs;
				a2        = 1 - alpha;
				break;
			case 'AllpassFO':
				var freq = conf.freq;
				var omega = 2 * Math.PI * freq / fs;
				var tn = Math.tan(omega / 2);
				var alpha = (tn + 1) / (tn - 1);
				b0 = 1;
				b1 = alpha;
				b2 = 0;
				a0 = alpha;
				a1 = 1;
				a2 = 0;
				break;
			case 'Bandpass':
				var freq  = conf.freq;
				var omega = 2 * Math.PI * freq / fs;
				var sn    = Math.sin( omega );
				var cs    = Math.cos( omega );
				var alpha = 'q' in conf ? sn / ( 2 * conf.q ) : sn * Math.sinh( Math.log( 2 ) / 2 * conf.bandwidth * omega / sn );
				b0        = alpha;
				b1        = 0;
				b2        = -alpha;
				a0        = 1 + alpha;
				a1        = -2 * cs;
				a2        = 1 - alpha;
				break;
			case 'Free':
				a0 = 1;
				a1 = conf.a1;
				a2 = conf.a2;
				b0 = conf.b0;
				b1 = conf.b1;
				b2 = conf.b2;
				break;
			case 'GeneralNotch':
				var f_p   = conf.freq_pole;
				var f_z   = conf.freq_zero;
				var q_p   = conf.q_pole;
				var tn_z  = Math.tan( Math.PI * f_z / fs );
				var tn_p  = Math.tan( Math.PI * f_p / fs );
				var alpha = tn_p / q_p;
				var tn2_p = tn_p ** 2;
				var tn2_z = tn_z ** 2;
				var gain  = conf.normalize_at_dc ? tn2_p / tn2_z : 1;
				b0        = gain * ( 1 + tn2_z ) ;
				b1        = -2 * gain * ( 1 - tn2_z );
				b2        = gain * ( 1 + tn2_z );
				a0        = 1 + alpha + tn2_p;
				a1        = -2 + 2 * tn2_p;
				a2        = 1 - alpha + tn2_p;
				break;
			case 'Highpass':
				var freqHP  = conf.freq;
				var qHP     = conf.q;
				var omegaHP = 2 * Math.PI * freqHP / fs;
				var snHP    = Math.sin( omegaHP );
				var csHP    = Math.cos( omegaHP );
				var alphaHP = snHP / ( 2 * qHP );
				b0          = ( 1 + csHP ) / 2;
				b1          = -( 1 + csHP );
				b2          = ( 1 + csHP ) / 2;
				a0          = 1 + alphaHP;
				a1          = -2 * csHP;
				a2          = 1 - alphaHP;
				break;
			case 'HighpassFO':
				var freq  = conf.freq;
				var omega = 2 * Math.PI * freq / fs;
				var k     = Math.tan( omega / 2 );
				var alpha = 1 + k;
				a0        = 1;
				a1        = -( ( 1 - k ) / alpha );
				a2        = 0;
				b0        = 1 / alpha;
				b1        = -1 / alpha;
				b2        = 0;
				break;
			case 'Highshelf':
				var freqHS  = conf.freq;
				var gainHS  = conf.gain;
				var omegaHS = 2 * Math.PI * freqHS / fs;
				var amplHS  = Math.pow( 10, gainHS / 40 );
				var snHS    = Math.sin( omegaHS );
				var csHS    = Math.cos( omegaHS );
				var alphaHS, betaHS;
				if ( 'slope' in conf ) {
					var slopeHS = conf.slope;
					alphaHS     = snHS / 2 * Math.sqrt( ( amplHS + 1 / amplHS) * ( 1 / ( slopeHS / 12 ) - 1 ) + 2 );
					betaHS      = 2 * Math.sqrt( amplHS ) * alphaHS;
				} else {
					var qHS     = conf.q;
					betaHS      = snHS * Math.sqrt( amplHS ) / qHS;
				}
				b0          = amplHS * ( ( amplHS + 1 ) + ( amplHS - 1 ) * csHS + betaHS );
				b1          = -2 * amplHS * ( ( amplHS - 1 ) + ( amplHS + 1 ) * csHS );
				b2          = amplHS * ( ( amplHS + 1 ) + ( amplHS - 1 ) * csHS - betaHS );
				a0          = ( amplHS + 1 ) - ( amplHS - 1 ) * csHS + betaHS;
				a1          = 2 * ( ( amplHS - 1 ) - ( amplHS + 1 ) * csHS );
				a2          = ( amplHS + 1 ) - ( amplHS - 1 ) * csHS - betaHS;
				break;
			case 'HighshelfFO':
				var freqHSFO  = conf.freq;
				var gainHSFO  = conf.gain;
				var omegaHSFO = 2 * Math.PI * freqHSFO / fs;
				var amplHSFO  = Math.pow( 10, gainHSFO / 40 );
				var tnHSFO    = Math.tan( omegaHSFO / 2 );
				b0            = amplHSFO * tnHSFO + Math.pow( amplHSFO, 2 );
				b1            = amplHSFO * tnHSFO - Math.pow( amplHSFO, 2 );
				b2            = 0;
				a0            = amplHSFO * tnHSFO + 1;
				a1            = amplHSFO * tnHSFO - 1;
				a2            = 0;
				break;
			case 'LinkwitzTransform':
				var f0  = conf.freq_act;
				var q0  = conf.q_act;
				var qt  = conf.q_target;
				var ft  = conf.freq_target;
				var d0i = Math.pow( 2 * Math.PI * f0, 2 );
				var d1i = ( 2 * Math.PI * f0 ) / q0;
				var c0i = Math.pow( 2 * Math.PI * ft, 2 );
				var c1i = ( 2 * Math.PI * ft ) / qt;
				var fc  = ( ft + f0 ) / 2;
				var gn  = 2 * Math.PI * fc / Math.tan( Math.PI * fc / fs );
				var cci = c0i + gn * c1i + Math.pow( gn, 2 );
				b0      = ( d0i + gn * d1i + Math.pow( gn, 2 ) ) / cci;
				b1      = 2 * ( d0i - Math.pow( gn, 2 ) ) / cci;
				b2      = ( d0i - gn * d1i + Math.pow( gn, 2 ) ) / cci;
				a0      = 1;
				a1      = 2 * ( c0i - Math.pow( gn, 2 ) ) / cci;
				a2      = ( c0i - gn * c1i + Math.pow( gn, 2 ) ) / cci;
				break;
			case 'Lowpass':
				var freqLP  = conf.freq;
				var qLP     = conf.q;
				var omegaLP = 2 * Math.PI * freqLP / fs;
				var snLP    = Math.sin( omegaLP );
				var csLP    = Math.cos( omegaLP );
				var alphaLP = snLP / ( 2 * qLP );
				b0          = ( 1 - csLP ) / 2;
				b1          = 1 - csLP;
				b2          = ( 1 - csLP ) / 2;
				a0          = 1 + alphaLP;
				a1          = -2 * csLP;
				a2          = 1 - alphaLP;
				break;
			case 'LowpassFO':
				var freq  = conf.freq;
				var omega = 2 * Math.PI * freq / fs;
				var k     = Math.tan( omega / 2 );
				var alpha = 1 + k;
				a0        = 1;
				a1        = -( ( 1 - k ) / alpha );
				a2        = 0;
				b0        = k / alpha;
				b1        = k / alpha;
				b2        = 0;
				break;
			case 'Lowshelf':
				var freqLS  = conf.freq;
				var gainLS  = conf.gain;
				var omegaLS = 2 * Math.PI * freqLS / fs;
				var amplLS  = Math.pow( 10, gainLS / 40 );
				var snLS    = Math.sin( omegaLS );
				var csLS    = Math.cos( omegaLS );
				var alphaLS, betaLS;
				if ( 'slope' in conf ) {
					var slopeLS = conf.slope;
					alphaLS = snLS / 2 * Math.sqrt( ( amplLS + 1 / amplLS ) * ( 1 / ( slopeLS / 12 ) - 1 ) + 2 );
					betaLS  = 2 * Math.sqrt( amplLS ) * alphaLS;
				} else {
					var qLS = conf.q;
					betaLS  = snLS * Math.sqrt( amplLS ) / qLS;
				}
				b0 = amplLS * ( ( amplLS + 1 ) - ( amplLS - 1 ) * csLS + betaLS );
				b1 = 2 * amplLS * ( ( amplLS - 1 ) - ( amplLS + 1 ) * csLS );
				b2 = amplLS * ( ( amplLS + 1 ) - ( amplLS - 1 ) * csLS - betaLS );
				a0 = ( amplLS + 1 ) + ( amplLS - 1 ) * csLS + betaLS;
				a1 = -2 * ( ( amplLS - 1 ) + ( amplLS + 1 ) * csLS );
				a2 = ( amplLS + 1 ) + ( amplLS - 1 ) * csLS - betaLS;
				break;
			case 'LowshelfFO':
				var freqLSFO  = conf.freq;
				var gainLSFO  = conf.gain;
				var omegaLSFO = 2 * Math.PI * freqLSFO / fs;
				var amplLSFO  = Math.pow( 10, gainLSFO / 40 );
				var tnLSFO    = Math.tan( omegaLSFO / 2 );
				b0            = Math.pow( amplLSFO, 2 ) * tnLSFO + amplLSFO;
				b1            = Math.pow( amplLSFO, 2 ) * tnLSFO - amplLSFO;
				b2            = 0;
				a0            = tnLSFO + amplLSFO;
				a1            = tnLSFO - amplLSFO;
				a2            = 0;
				break;
			case 'Notch':
				var freq  = conf.freq;
				var omega = 2 * Math.PI * freq / fs;
				var sn    = Math.sin( omega );
				var cs    = Math.cos( omega );
				var alpha = 'q' in conf ? sn / ( 2 * conf.q ) : sn * Math.sinh( Math.log( 2 ) / 2 * conf.bandwidth * omega / sn );
				b0        = 1;
				b1        = -2 * cs;
				b2        = 1;
				a0        = 1 + alpha;
				a1        = -2 * cs;
				a2        = 1 - alpha;
				break;
			case 'Peaking':
				var freqPK  = conf.freq;
				var gainPK  = conf.gain;
				var omegaPK = 2 * Math.PI * freqPK / fs;
				var snPK    = Math.sin( omegaPK );
				var csPK    = Math.cos( omegaPK );
				var amplPK  = Math.pow( 10, gainPK / 40 );
				var alphaPK;
				if ( 'q' in conf ) {
					var qPK = conf.q;
					alphaPK = snPK / ( 2 * qPK );
				} else {
					var bandwidthPK = conf.bandwidth;
					alphaPK         = snPK * Math.sinh( Math.log( 2 ) / 2 * bandwidthPK * omegaPK / snPK );
				}
				b0 = 1 + ( alphaPK * amplPK );
				b1 = -2 * csPK;
				b2 = 1 - ( alphaPK * amplPK );
				a0 = 1 + ( alphaPK / amplPK );
				a1 = -2 * csPK;
				a2 = 1 - ( alphaPK / amplPK );
				break;
		}
		this.fs = fs;
		this.a1 = a1 / a0;
		this.a2 = a2 / a0;
		this.b0 = b0 / a0;
		this.b1 = b1 / a0;
		this.b2 = b2 / a0;
	}

	complexGain( freq, removeDelay = false ) {
		var fL = freq.length;
		var A  = [];
		var a, b, z, z_1, z_2;
		for ( var i = 0; i < fL; i++ ) {
			z   = Complex( 0, 1 )
							.mul( 2 * Math.PI * freq[ i ] / this.fs )
							.exp();
			z_1 = z.pow( -1 );
			z_2 = z.pow( -2 );
			b   = z_1
					.mul( this.b1 )
					.add( this.b0 )
					.add( z_2.mul( this.b2 ) );
			a   = z_1
					.mul( this.a1 )
					.add( 1 )
					.add( z_2.mul( this.a2 ) );
			A.push( b.div( a ) );
		}
		return [ freq, A ];
	}
}
class BiquadCombo extends BaseFilter {
	ButterwQ( order ) {
		var odd     = order % 2 > 0;
		var n_so    = Math.floor( order / 2 );
		var qvalues = [];
		for ( var n = 0; n < n_so; n++ ) {
			var q = 1 / ( 2 * Math.sin( ( Math.PI / order ) * ( n + 1 / 2 ) ) );
			qvalues.push( q );
		}
		if ( odd ) {
			qvalues.push( -1 );
		}
		return qvalues;
	}

	constructor(conf, fs) {
		super();
		this.ftype = conf.type;
		switch ( this.ftype ) {
			case 'ButterworthHighpass':
			case 'ButterworthLowpass':
			case 'LinkwitzRileyHighpass':
			case 'LinkwitzRileyLowpass':
				this.order  = conf.order;
				this.freq   = conf.freq;
				this.fs     = fs;
				this.typeSO = this.ftype.includes( 'Hi' ) ? 'Highpass' : 'Lowpass';
				this.typeFO = this.typeSO +'FO';
				var qvalues;
				switch ( this.ftype ) {
					case 'LinkwitzRileyHighpass':
						var q_temp = this.ButterwQ( this.order / 2 );
						qvalues = ( this.order / 2 ) % 2 > 0 ? q_temp.slice( 0, -1 ).concat( q_temp ).concat( [ 0.5 ] ) : q_temp.concat( q_temp );
						break;
					case 'LinkwitzRileyLowpass':
						var q_temp = this.ButterwQ(this.order / 2);
						qvalues = ( this.order / 2 ) % 2 > 0 ? q_temp.slice( 0, -1 ).concat( q_temp ).concat( [ 0.5 ] ) : q_temp.concat( q_temp );
						break;
					case 'ButterworthHighpass':
						qvalues = this.ButterwQ(this.order);
						break;
					case 'ButterworthLowpass':
						qvalues = this.ButterwQ( this.order );
						break;
				}
				this.biquads = qvalues.map( q => new Biquad( q >= 0 ? { freq: this.freq, q: q, type: this.typeSO } : { freq: this.freq, type: this.typeFO }, this.fs ) );
				break;
			case 'FivePointPeq':
				this.biquads = [
					new Biquad( { freq: conf.fls, q: conf.qls, gain: conf.gls, type: 'Lowshelf' },  fs ),
					new Biquad( { freq: conf.fp1, q: conf.qp1, gain: conf.gp1, type: 'Peaking' },   fs ),
					new Biquad( { freq: conf.fp2, q: conf.qp2, gain: conf.gp2, type: 'Peaking' },   fs ),
					new Biquad( { freq: conf.fp3, q: conf.qp3, gain: conf.gp3, type: 'Peaking' },   fs ),
					new Biquad( { freq: conf.fhs, q: conf.qhs, gain: conf.ghs, type: 'Highshelf' }, fs )
				];
				break;
			case 'GraphicEqualizer':
				var bands     = conf.gains.length;
				var f_min     = conf.freq_min || 20;
				var f_max     = conf.freq_max || 20000;
				var f_min_log = Math.log2( f_min );
				var f_max_log = Math.log2( f_max );
				var bw        = ( f_max_log - f_min_log ) / bands;
				this.biquads  = conf.gains.map( ( gain, band ) => {
					var freq_log = f_min_log + ( band + 0.5 ) * bw;
					var freq     = Math.pow( 2, freq_log );
					return new Biquad( { freq: freq, bandwidth: bw, gain: gain, type: 'Peaking' }, fs );
				});
				break;
			case 'Tilt':
				var gain_low = -conf.gain / 2;
				var gain_high = conf.gain / 2;
				this.biquads = [
					new Biquad( { freq:  110, q: 0.35, gain: gain_low,  type: 'Lowshelf' },  fs ),
					new Biquad( { freq: 3500, q: 0.35, gain: gain_high, type: 'Highshelf' }, fs )
				];
				break;
		}
	}

	complexGain( freq, removeDelay = false ) {
		var A = new Array( freq.length ).fill( 1 );
		this.biquads.forEach( bq => {
			var [ _, Atemp ] = bq.complexGain( freq );
			A                = A.map( ( a, i ) => Atemp[ i ].mul( a ) );
		});
		return [ freq, A ];
	}
}
class Conv {
	constructor( conf, fs ) {
		if ( ! conf ) conf = { values: [1] };
		this.conf = conf;
		this.fs   = fs;
	}
	
	findPeak() {
		var maxVal = Math.abs( this.impulse[ 0 ] );
		var idx    = 0;
		for ( var i = 1; i < this.impulse.length; i++ ) {
			var absVal = Math.abs( this.impulse[ i ] );
			if ( absVal > maxVal ) {
				maxVal = absVal;
				idx    = i;
			}
		}
		return idx;
	}

	complexGain( freq, removeDelay = false ) {
		var impulselen = this.impulse.length;
		var npoints    = Math.pow( 2, Math.ceil( Math.log2( impulselen ) ) );
		if ( npoints < 1024 ) npoints = 1024;
		var impulse    = this.impulse.slice();
		var padding    = new Array( npoints - impulselen ).fill( 0 );
		impulse        = impulse.concat( padding );
		var f_fft      = Array.from( { length: npoints / 2 }, ( _, n ) => this.fs * n / npoints );
		var cut        = this.fft_im.slice( 0, npoints / 2 );
		if ( removeDelay ) {
			var maxidx = this.findPeak();
			var cutL   = cut.length;
			var Avec   = [];
			for ( i = 0; i < cutL; i++ ) {
				Avec.push( Complex( 0, 1 )
								.neg()
								.mul( Math.PI * i * maxidx / ( npoints / 2 ) )
								.exp()
								.mul( cut[ i ] ) );
			}
		}
		if ( freq !== null ) {
			var interpolated = this.interpolatePolar( Avec, f_fft, freq );
			return [ freq, interpolated ];
		}
		return [ f_fft, removeDelay ? Avec : cut ];
	}

	interpolate( y, xold, xnew ) {
		var ynew = [];
		for ( var x of xnew ) {
			var idx = ( y.length * x ) / xold[ xold.length - 1 ];
			var i1  = Math.floor( idx );
			var i2  = i1 + 1;
			if ( i1 >= y.length ) {
				i1 = y.length - 1;
			}
			if ( i2 >= y.length ) {
				i2 = i1;
			}
			var fract  = idx - i1;
			var newval = ( 1 - fract ) * y[ i1 ] + fract * y[ i2 ];
			ynew.push( newval );
		}
		return ynew;
	}

	interpolatePolar( y, xold, xnew ) {
		var y_magn = y.map( yval => yval.abs() );
		var y_ang  = mapPhase( y );
		y_ang      = unwrapPhase( y_ang, 270 ).map( yval => Math.PI * yval / 180 );
		var y_magn_interp = this.interpolate (y_magn, xold, xnew );
		var y_ang_interp  = this.interpolate( y_ang, xold, xnew );
		return y_magn_interp.map( ( r, i ) => Complex( { r, phi: y_ang_interp[ i ] } ) );
	}

	gainAndPhase( freq, removeDelay = false, callback ) {
		bash( "camilla_audiofileread.py '"+ JSON.stringify( this.conf ) +"'", data => {
			this.impulse        = data[ 0 ];
			this.fft_im         = [];
			data[ 1 ].forEach( ri => this.fft_im.push( Complex( ri ) ) );
			var [ f_fft, Avec ] = this.complexGain(null, removeDelay);
			Avec                = this.interpolatePolar(Avec, f_fft, freq);
			var ga              = mapGain( Avec );
			var ph              = mapPhase( Avec );
			var [ fg, gr ]      = calcGroupDelay( freq, ph );
			callback( fg, ga, gr, ph, this.impulse );
		}, 'json' );
	}

	getImpulse() {
		var t = Array.from( { length: this.impulse.length }, ( _, n ) => n / this.fs );
		return [ t, this.impulse ];
	}
}
class Delay extends BaseFilter {
	constructor( conf, fs ) {
		super();
		this.fs = fs;
		var unit = conf.unit || 'ms';
		if ( unit === 'ms' ) {
			this.delaySamples = conf.delay / 1000 * fs;
		} else if ( unit === 'mm' ) {
			this.delaySamples = conf.delay / 1000 * fs / 343;
		} else if ( unit === 'samples' ) {
			this.delaySamples = conf.delay;
		} else {
			throw new Error( 'Unknown unit: '+ unit );
		}

		this.subsample = conf.subsample === true;
		if ( this.subsample ) {
			this.delayFullSamples = Math.floor( this.delaySamples );
			this.fraction = this.delaySamples - this.delayFullSamples;
			this.a1 = 1 - this.fraction;
			this.a2 = 0;
			this.b0 = 1 - this.fraction;
			this.b1 = 1;
			this.b2 = 0;
		} else {
			this.delayFullSamples = Math.round( this.delaySamples );
		}
	}

	complexGain( freq, removeDelay = false ) {
		var zvec = freq.map( f => {
			Complex( 0, 1 )
					.mul( 2 * Math.PI * f / this.fs )
					.exp();
		} );
		var A;
		if ( this.subsample ) {
			var a, b, z, z_1, z_2;
			A = zvec.map( z => {
				z   = Complex( 0, 1 )
								.mul( 2 * Math.PI * freq[ i ] / this.fs )
								.exp();
				z_1 = z.pow( -1 );
				z_2 = z.pow( -2 );
				b   = z_1
						.mul( this.b1 )
						.add( this.b0 )
						.add( z_2.mul( this.b2 ) );
				a   = z_1
						.mul( this.a1 )
						.add( 1 )
						.add( z_2.mul( this.a2 ) );
				b.div( a );
			} );
		} else {
			A = new Array( zvec.length ).fill( 1 );
		}
		if ( ! removeDelay ) {
			var delayS = this.delayFullSamples / this.fs;
			A = A.map( ( val, i ) => {
				Complex( 0, 1 )
						.mul( 2 * Math.PI * freq[ i ] * delayS )
						.exp()
						.mul( val );
			} );
		}
		return [ freq, A ];
	}
}
class DiffEq extends BaseFilter {
	constructor( conf, fs ) {
		super();
		this.fs = fs;
		this.a = conf.a;
		this.b = conf.b;
		if ( this.a.length === 0 ) this.a = [1];
		if ( this.b.length === 0 ) this.b = [1];
	}

	complex_gain( freq, remove_delay = false ) {
		var zvec = freq.map( f => {
			Complex( 0, 1 )
					.mul( 2 * Math.PI * f / this.fs)
					.exp();
		} );
		var A1   = new Array( freq.length ).fill( 0 );
		this.b.forEach( ( bn, n ) => {
			A1 = A1.map( ( a1, i ) => {
				zvec[ i ]
					.pow( -n )
					.mul( bn )
					.add( a1 );
			} );
		} );
		var A2   = new Array( freq.length ).fill( 0 );
		this.a.forEach( ( an, n ) => {
			A2 = A2.map( ( a2, i ) => {
				zvec[ i ]
					.pow( -n )
					.mul( an )
					.add( a2 );
			} );
		} );
		var A    = A1.map( ( a1, i ) => a1.div( A2[ i ] ) );
		return [ freq, A ];
	}
}
class Gain extends BaseFilter {
	constructor( conf ) {
		super();
		this.gain     = conf.gain;
		this.inverted = conf.inverted === true;
		this.scale    = conf.scale || 'dB';
	}

	complexGain( freq, removeDelay = false ) {
		var sign = this.inverted ? -1 : 1;
		var gain;
		if ( this.scale === 'dB' ) {
			gain = Math.pow( 10, this.gain / 20 ) * sign;
		} else {
			gain = this.gain * sign;
		}
		var A = new Array( freq.length ).fill( gain );
		return [ freq, A ];
	}
}
class Loudness extends BaseFilter {
	constructor( conf, fs, volume ) {
		super();
		var relVol    = volume - conf.reference_level;
		var relBoost  = -relVol / 20;
		relBoost      = Math.max( 0, Math.min( 1, relBoost ) );

		var highBoost = relBoost * conf.high_boost;
		var lowBoost  = relBoost * conf.low_boost;

		if ( conf.attenuate_mid ) {
			var maxGain  = Math.max( highBoost, lowBoost );
			this.midGain = Math.pow( 10, -maxGain / 20 );
		} else {
			this.midGain = 1;
		}

		var lsConf    = new Biquad(
			{ freq: 70, slope: 12, gain: lowBoost, type: 'Lowshelf' }, 
			fs
		);
		var hsConf    = new Biquad(
			{ freq: 3500, slope: 12, gain: highBoost, type: 'Highshelf' }, 
			fs
		);
		this.biquads  = [ lsConf, hsConf ];
	}

	complexGain( freq, removeDelay = false ) {
		var A = new Array( freq.length ).fill( this.midGain );
		this.biquads.forEach( bq => {
			var [ _, Atemp ] = bq.complexGain( freq );
			A                = A.map( ( a, i ) => Atemp[ i ].mul( a ) );
		});
		return [ freq, A ];
	}
}
