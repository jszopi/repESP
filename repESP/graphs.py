import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import field_comparison


def plot3d(esp_field, diff, dist, exclusion_dist=0, rand_skim=0.02):
    field_comparison._check_grids(esp_field, diff, dist)
    if esp_field.field_type != 'esp':
        raise TypeError("The field passed was of type '{0}' but 'esp' was "
                        "expected.".format(esp_field.field_type))
    if 'diff' not in diff.field_type:
        raise TypeError("The field passed was of type '{0}' but one of the "
                        "'diff' types was expected.".format(diff.field_type))
    if 'dist' not in dist.field_type:
        raise TypeError("The field passed was of type '{0}' but one of the "
                        "'dist' types was expected.".format(dist.field_type))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.set_xlabel('Distance ' + dist.field_type)
    ax.set_ylabel('ESP value')
    ax.set_zlabel('ESP ' + diff.field_type)

    dist, esp_field, diff = field_comparison.skim(rand_skim, dist, esp_field,
                                                  diff)
    dist, esp_field, diff = field_comparison.filter_by_dist(
        exclusion_dist, dist, esp_field, diff)
    # Flatten and remove NANs
    dist, esp_field, diff = map(field_comparison._flatten_no_nans, (
        dist, esp_field, diff))

    ax.scatter(dist, esp_field, diff)
    plt.show()
    plt.close()
