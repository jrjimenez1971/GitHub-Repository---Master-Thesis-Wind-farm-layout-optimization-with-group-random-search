# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 22:44:43 2025

@author: jrjim
"""

import py_wake
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path

from py_wake.examples.data.hornsrev1 import Hornsrev1Site
from py_wake.examples.data.hornsrev1 import wt_x, wt_y

def create_mesh_border(points):
    """
    Creates a matplotlib.path.Path object from a list of (x, y) points.

    Args:
        points: A list of tuples representing the border points.

    Returns:
        A matplotlib.path.Path object.
    """
    codes = [Path.MOVETO] + [Path.LINETO] * (len(points) - 1)
    return Path(points, codes)

def filter_points_inside_border(border_path, points):
    """
    Filters a list of (x, y) points to identify those inside or outside a given border.

    Args:
        border_path: A matplotlib.path.Path object representing the border.
        points: A list or NumPy array of (x, y) coordinates to filter.

    Returns:
        A tuple containing two lists: (inside_points, outside_points).
    """

    inside_points = []
    outside_points = []

    for x, y in points:
        if border_path.contains_point((x, y)):
            inside_points.append((x, y))
        else:
            outside_points.append((x, y))

    return inside_points, outside_points

def visualize_results(border_points, inside_points, outside_points):
    """
    Visualizes the border and the filtered points.

    Args:
        border_points: A list of (x, y) coordinates defining the border.
        inside_points: A list of (x, y) coordinates inside the border.
        outside_points: A list of (x, y) coordinates outside the border.
    """
    border_x, border_y = zip(*border_points)
    plt.plot(border_x, border_y, 'r-', label='Border')

    if inside_points:
        inside_x, inside_y = zip(*inside_points)
        plt.scatter(inside_x, inside_y, c='g', marker='o', label='Inside')

    if outside_points:
        outside_x, outside_y = zip(*outside_points)
        plt.scatter(outside_x, outside_y, c='b', marker='x', label='Outside')

    plt.legend()
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Points Inside/Outside Border')
    plt.grid(True)
    plt.axis('equal') #ensure consistent scales
    plt.show()

if __name__ == "__main__":
    # Define the border points (example: a simple polygon)
    # border_points = [(0, 0), (2, 4), (5, 3), (4, 0), (0, 0)]
    border_points = [((wt_x[0]-25), (wt_y[0]+25)), ((wt_x[7]-25), (wt_y[7]-25)), ((wt_x[79]+25), (wt_y[79]-25)), ((wt_x[72]+25, wt_y[72]+25)), ((wt_x[0]-25, wt_y[0]+25))]


    # Create the border path
    border_path = create_mesh_border(border_points)

    # Generate some random test points
    num_points = 80
    # test_points = np.random.rand(num_points, 2) * 6  # Random points within a range
    test_points = []
    for i in range (0,79):
        wt_point = wt_x[i],wt_y[i]    
        test_points.append(wt_point)

    # Filter the points
    inside, outside = filter_points_inside_border(border_path, test_points)

    # Visualize the results
    visualize_results(border_points, inside, outside)