Rscript -e "require(knitr) ; require(markdown) ;knit('spatial_resolution.Rmd', 'spatial_resolution.md'); markdownToHTML('spatial_resolution.md', 'spatial_resolution.html');"
rm spatial_resolution.md
