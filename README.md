# StructureFromMotion

 I developed a 3D reconstruction algorithm as part of a computer vision course.
The goal was to reconstruct a scene in 3D using two successive images from the KITTI Vision Benchmark.
The project was developed on Matlab with the Image Processing Toolbox and Peter Kovesi's toolbox (http://www.peterkovesi.com/matlabfns/index.html).


Algorithm :

- Interest points detection using the Harris detector on both images
- Matching by correlation of the detected points
- Robust estimation of the fundamental matrix
- Projection matrices reconstruction in the two images
- Computation of the 3D structure by triangulation


#### Example of reconstruction
![alt tag](images/reconstruction.png?raw=true "Reconstruction")


#### Images
![alt tag](images/image1.png?raw=true "Image 1")
![alt tag](images/image2.png?raw=true "Image 2")


