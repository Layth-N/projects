#!/usr/bin/env python
"""mapper.py"""

import sys
from math import sqrt

#get initial centroids from a txt file and add them in an array
def getCentroids():
    centroids = []
    with open("C:\KMeanClustering\Centroids5.txt", 'r') as file: #reading the txt file from its location
        for line in file:
            line = line.strip()
            centroid = list(map(float, line.split(','))) #splitting the values
            centroids.append(centroid) #creating centroids array
    return centroids

#create clusters based on initial centroids
def createClusters(point, centroids):
    min_distance = float('inf') #initialising minimum distance
    cluster_id = None
    for i, centroid in enumerate(centroids):
        distance = sqrt(sum((x - y) ** 2 for x, y in zip(point, centroid))) #calculating the distance of the current point with each centroid
        if distance < min_distance: #assigning the point to the centroid by comparing the distances
            min_distance = distance
            cluster_id = i
    return cluster_id, point

if __name__ == "__main__":
    centroids = getCentroids()
    for line in sys.stdin: #going through each line/point
        line = line.strip()
        try:  
            point = list(map(float, line.split(','))) #splitting the points and turning them into floats
        except ValueError:  
        #point was not a number/deformed, so silently  
        #ignore/discard this line  
            continue
        cluster_id, assigned_point = createClusters(point, centroids) #combining the cluster ID with the point
        print(f"{cluster_id}\t{','.join(map(str, assigned_point))}")
