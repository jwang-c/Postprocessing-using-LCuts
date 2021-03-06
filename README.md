# Postprocessing-using-LCuts
LCuts: Linear Clustering of Bacteria using Recursive Graph Cuts. This repo will showcase the post-processing procedure on under-segmented cluster using LCuts for 3D bacterial dataset.

## Related papers
Please cite the following papers if you are referring the code. And please let me know if you have any questions. Thank you!
Email: jiewang@virginia.edu

[1] Wang J, Batabyal T, Zhang M, Zhang J, Aziz A, Gahlmann A, Acton ST. Lcuts: Linear Clustering of Bacteria Using Recursive Graph Cuts. In 2019 IEEE International Conference on Image Processing (ICIP) 2019 Sep 22 (pp. 1575-1579). IEEE.

[2] Zhang M, Zhang J, Wang Y, Wang J, Achimovich AM, Acton ST, Gahlmann A. Non-Invasive Single-Cell Morphometry and Tracking in Living Bacterial Biofilms. Submitting soon.

For other codes related to [2], please see [BCM3D](https://github.com/GahlmannLab/BCM3D).

## Code description
There are four major steps included in the main function.
- Data filtering for selecting the region of interest.
- Point cloud data generation by constrained medial axis extraction.
- Linear clustering using LCuts.
- Biofilm reconstrcution with geometrical models.

## Run the code
Just run the [demo.m](demo.m) file. 

The [post_processing_using_lcuts.m](post_processing_using_lcuts.p) file introduces the whole algorithm. Other detailed instructions can be found in each function.

Currently, parameters are chosen based on the image and bacterial cell information. Please look out for later updates for further information.

A Matlab version of R2018a or higher might be required.

The details of the source code are currenly private, but they will be public soon!

## License
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments
Thanks for the support from Virginia Image and Video Analysis (VIVA) Lab, University of Virginia.

Thanks for the data and biological/chemical knowledge support from the Cell imaging at Nano scale Lab, University of Virginia.
