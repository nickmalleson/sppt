Rscript -e "require(knitr) ; require(markdown) ;knit('spatial_resolution2.Rmd', 'spatial_resolution2.md'); markdownToHTML('spatial_resolution2.md', 'spatial_resolution2.html');"
rm spatial_resolution2.md
