# Radar target generation and detection project
## Implementation steps for the 2D CFAR process.
Instead of using for loops to calculate the thresholds, I use convolution to calculate those values. <br>
The threshold of each cell is the average noise value from the training
cells. Use a convolution with a customized kernel to produce the result.<br>
Here is an example of kernel which Tr=Gr=Td=Gd=1<br>

1/16, 1/16, 1/16, 1/16, 1/16,<br>
1/16, 0,    0,    0,    1/16,<br>
1/16, 0,    0,    0,    1/16,<br>
1/16, 0,    0,    0,    1/16,<br>
1/16, 1/16, 1/16, 1/16, 1/16,<br>

## Selection of Training, Guard cells and offset.
The numbers of training cells for range and doppler are selected as 10 and 8 to produce a smooth threshold.<br>
The numbers of guard cells are both equal to 4.<br>
The offset value is adjusted to 8 to produce similar results like the project walkthrough.

## Steps taken to suppress the non-thresholded cells at the edges.
The process above will generate a thresholded block, which is smaller than the Range Doppler Map as the CUT cannot be located at the edges of the matrix. Hence, few cells will not be thresholded. To keep the map size the same, set those values to 0.
