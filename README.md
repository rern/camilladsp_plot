Convert from [pycamilladsp-plot](https://github.com/HEnquist/pycamilladsp-plot/blob/master/camilladsp_plot/)
- `camilladsp_plot.js` << `filters.py`
	- require `complex.js` - math operations >> serial of complex function
  	```js
    Complex( 0, 1 )
  			.mul( 2 * Math.PI * f / this.fs )
  			.exp()
    ```
  - add `fft()` for `Conv` `complexGain()`
- `camilla.js`:
	- `graph.filter()`, `graph.pipeline()` << `eval_filterconfig.py`
- `camilla_audiofilread.py` << `audiofilread.py`
	- modify for `graph.filter()` - `Conv`
- `camilla_cooley_tukey.py` < `cooley_tukey.py`
  - (as is)
