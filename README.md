# overlapping-circles
Python code that calculates the total area of any number of overlapping (or non-overlapping) circles.  
The list `circles` should be a 2D list of 2x1 numpy arrays (center coordinates) paired with each circle's radius.  

For example, the pre-inputted circles in `main.py` generate the following:  
![image](https://user-images.githubusercontent.com/107213996/212503842-14e63d13-1f0b-42e9-b418-abfee27b873d.png)
  
The algorithm I have used is from:  
https://stackoverflow.com/a/24042372  

Behind the scenes, the area above is partitioned into polygons and circle sectors:  
![image](https://user-images.githubusercontent.com/107213996/212503877-bbe41dcc-7c49-43cd-a081-44818c281a78.png)  
So that no area is overlapping.  
  
We then perform some python magic and find the areas of each subsection, and add them together. 
