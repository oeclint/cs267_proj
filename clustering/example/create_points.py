import numpy as np
import os
import argparse
#import matplotlib.pyplot as plt

def create_points(error, num_samples, filename):
    if os.path.isfile(filename):
        os.remove(filename)

#    for i in [10,20]:
#        rad = np.random.randn(num_samples // 2) * 2 * np.pi
#        points_x, points_y = i*np.cos(rad), i*np.sin(rad)
#        points_x += np.random.uniform(-1*error, error, size = num_samples//2)
#        points_y += np.random.uniform(-1*error, error, size = num_samples//2)

    for i in [-10, 10]:
        rad = np.random.randn(num_samples // 2) * 2 * np.pi
        points_x, points_y = np.cos(rad) + i, np.sin(rad) + i
        with open(filename, 'a+') as f:
            for j in range(num_samples // 2):
                f.write(str(round(points_x[j], 4)) + '\t' + str(round(points_y[j], 4)) + '\t' + str((i/10)-1) + '\n')
#        plt.scatter(points_x, points_y)
#        plt.show()

def main():
    parser = argparse.ArgumentParser(description = "Create points")
    parser.add_argument('num_samples', type = int)
    parser.add_argument('error', type = float)
    parser.add_argument('filename', type = str)
    args = parser.parse_args()
    error = args.error
    num_samples = args.num_samples
    filename = args.filename
    #print(error, num_samples)
    create_points(error, num_samples, filename)

if __name__ == '__main__':
    main()
