This branch contains the input files used for the DynoPath experiemnt runs.
File names define the specific grid the files relate to.
Example: FileName = 42X122_In20S_Spread_Out60_Spred.txt
42X122 = Grid size of 40X120. Related as 42X122 for extra 2 rows and 2 columns for input and output cells excluded from the grid size.
In20S_Spread = 20 Input cells spread on the South side of the grid
Out60_Spred = 60 Outpute cells spread along the other three edges of the grid

The Data in the file contains 4 rows:
42,122 = X,Y dimensions
(1,0),(3,0),(5,0).... = locations of the input cells in (x,y) format where (0,0) is the left bottom corner
(2,121,0),(6,121,1),(10,121,2)... = locations of the output cells in (x,y,z) format where x,y define the cell location and z defines the item type (can be more the output cell per item type)
(0,0),(0,1),(0,2)... = locations of inactive cells
