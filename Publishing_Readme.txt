For publishing the documentation, in general, use

 opts.showCode = false
 opts.evalCode = false
 opts.format = 'html'
 publish('Scyllarus.m',opts);

For the function_list.m and example_code.m use

 publish('function_list.m');

Note that the function_list.m file is to be updated every time a function is added. 
The URLs on the file need to be made consistent with any website the documentation is on. 
Also, remember to clean up the directory structure so as to keep it kneat.