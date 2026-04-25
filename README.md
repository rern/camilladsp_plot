### For rAudio
Converted from [pycamilladsp-plot](https://github.com/HEnquist/pycamilladsp-plot/blob/master/camilladsp_plot/)
- `/srv/http/assets/js/camilla.js` << `eval_filterconfig.py`
	- `graph.filter()`
	- `graph.pipeline()`
- `/srv/http/assets/js/plugin/camilladsp_plot.js` << `filters.py`
	- add `fft()` for `Conv` - `complexGain()`
- `/srv/http/assets/js/plugin/complex.js` - required for math operations
  	```js
	// serial of complex function
    Complex( 0, 1 )
  			.mul( 2 * Math.PI * f / this.fs )
  			.exp()
    ```
- `/srv/http/bash/settings/camilla_audiofilread.py` << `audiofilread.py`
	- modify for `Conv` - `graph.filter()`
- `/srv/http/bash/settings/camilla_cooley_tukey.py` < `cooley_tukey.py`
	- copy
