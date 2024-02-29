#/bin/bash

cd PDF

pdfseparate allplot_$(date -d "-1 days" +%Y-%m-%d).pdf merged-%02d-1.pdf
pdfseparate allplot_$(date +%Y-%m-%d).pdf              merged-%02d-2.pdf

pdfjam --landscape --nup 2x1 *merged*

mv merged*pdfjam.pdf diffplot_$(date +%Y-%m-%d).pdf

\rm -f *merged*

cd -
 



