#!/usr/bin/env python
"""reducer.py"""


import sys

with open("C:\KMeanClustering\Centroids5.txt", 'w') as file: #emptying the centroids file
    pass

def calculateNewCentroids():
    current_centroid = None
    current_sum = None  # sum of points for each cluster
    cluster_count = 0  # count of points for each cluster
    
    with open("C:\KMeanClustering\Centroids5.txt", 'a') as file:
    # input comes from STDIN
        for line in sys.stdin:
            line = line.strip()
            cluster_id, point = line.split('\t')
            cluster_id_int = int(cluster_id)
            point_split = list(map(float, point.strip().split(',')))

            if current_centroid == cluster_id_int:
                current_sum = [x + y for x, y in zip(current_sum, point_split)]
                cluster_count += 1
            else:
                if current_centroid is not None:
                    new_centroid = [x / cluster_count for x in current_sum]
                    centroid_string = ','.join(map(str, new_centroid))
                    print(centroid_string)
                    file.write(centroid_string + '\n')
                
                current_centroid = cluster_id_int
                current_sum = point_split
                cluster_count = 1
        if current_centroid is not None:
                    new_centroid = [x / cluster_count for x in current_sum]
                    centroid_string = ','.join(map(str, new_centroid))
                    print(centroid_string)
                    file.write(centroid_string)

if __name__ == "__main__":
    calculateNewCentroids()
