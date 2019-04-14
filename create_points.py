import numpy as np
import os
import argparse

def create_points(error, num_samples):
    if os.path.isfile('points.txt'):
        os.remove('points.txt')

    for i in [1,2]:
        rad = np.random.randn(num_samples // 2) * 2 * np.pi
        points_x, points_y = i*np.cos(rad), i*np.sin(rad)
        points_x += np.random.uniform(-1*i*error, i*error, size = num_samples//2)
        points_y += np.random.uniform(-1*i*error, i*error, size = num_samples//2)
        with open('points.txt', 'a+') as f:
            for j in range(num_samples // 2):
                f.write(str(round(points_x[j], 4)) + '\t' + str(round(points_y[j], 4)) + '\t' + str(i-1) + '\n')
        #plt.scatter(points_x, points_y)
        #plt.show()

def main():
    parser = argparse.ArgumentParser(description = "Create points")
    parser.add_argument('num_samples', type = int)
    parser.add_argument('error', type = float)
    args = parser.parse_args()
    error = args.error
    num_samples = args.num_samples
    #print(error, num_samples)
    create_points(error, num_samples)

if __name__ == '__main__':
    main()
