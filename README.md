Convert from [pycamilladsp-plot](https://github.com/HEnquist/pycamilladsp-plot/blob/master/camilladsp_plot/)
- `/srv/http/assets/js/plugin/camilladsp_plot.js` << `filters.py`
	- require `complex.js` - math operations
  	```js
	// serial of complex function
    Complex( 0, 1 )
  			.mul( 2 * Math.PI * f / this.fs )
  			.exp()
    ```
  - add `fft()` for `Conv` `complexGain()`
- `/srv/http/assets/js/camilla.js`:
	- `graph.filter()`, `graph.pipeline()` << `eval_filterconfig.py`
- `/srv/http/bash/settings/camilla_audiofilread.py` << `audiofilread.py`
	- modify for `graph.filter()` - `Conv`
- `/srv/http/bash/settings/camilla_cooley_tukey.py` < `cooley_tukey.py`
	- copy
