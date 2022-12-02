import helper as my
import operator
import math
import statistics

# =============================================================================
# SECTION 1
# =============================================================================
def distance(v, w):
    """
    Calculates the distance between m-dimensional vectors v and w, where
    the distance is the square root of the sum of squared differences for each
    point i in m.

    Parameters
    ----------
    v : LIST
        Vector V for i = 1, 2, ..., m.
    w : LIST
        Vector V for i = 1, 2, ..., m.

    Returns
    -------
    FLOAT
        Square root of sum of squared difference between v_i and w_i.

    """
    diff = [i**2 for i in list(map(operator.sub, v, w))]
    return math.sqrt(sum(diff))


def distances(v, centres):
    """
    For data point v, distances() calculates the distance between v and all of
    centres.

    Parameters
    ----------
    v : LIST
        m-dimensional vector.
    centres : LIST of LISTS
        List of centres of length n, each containing a centre j of length m.

    Returns
    -------
    FLOAT
        The minimum distance between v and any centre.
    centre : LIST
        Returns the centre for which the distance is minimised.

    """
    d = [distance(v, w) for w in centres]
    centre = d.index(min(d))
    return min(d), centre


def max_distance(dataset, centres):
    """
    Calculates the farthest distance from any datapoint and its assigned
    centre.

    Parameters
    ----------
    dataset : LIST of LISTS
        Dataset of length n containing n vectors of length m.
    centres : LIST of LISTS
        List of k centres of length m.

    Returns
    -------
    FLOAT
        Maximum distance.
    index : INT
        Index of datapoint which has max distance.
    centre : TYPE
        Datapoint which has max distance.

    """
    d = []
    for data in dataset:
        distances = [distance(data, centre) for centre in centres]
        d.append(min(distances))
    index = d.index(max(d))
    return max(d), index, dataset[index]


def farthest_first_traversal(data, k):
    """
    

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    k : TYPE
        DESCRIPTION.

    Returns
    -------
    centres : TYPE
        DESCRIPTION.

    """
    centres = [data[0]]
    while len(centres) < k:
        temp = max_distance(data, centres)
        centres.append(data[temp[1]])
    return centres


def distortion(dataset, centres):
    """
    

    Parameters
    ----------
    dataset : TYPE
        DESCRIPTION.
    centres : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    d = []
    for data in dataset:
        d.append(distances(data, centres)[0]**2)
    return round(sum(d)/len(dataset), 3)


def gravity(cluster):
    """
    

    Parameters
    ----------
    cluster : TYPE
        DESCRIPTION.

    Returns
    -------
    res : TYPE
        DESCRIPTION.

    """
    # n = len(cluster)
    # res = [sum(x)/n for x in zip(*cluster)]
    res = [round(statistics.mean(x), 3) for x in zip(*cluster)]
    return res


def centres_to_clusters(dataset, centres):
    """
    

    Parameters
    ----------
    dataset : TYPE
        DESCRIPTION.
    centres : TYPE
        DESCRIPTION.

    Returns
    -------
    clusters : TYPE
        DESCRIPTION.

    """
    clusters = {new_list: [] for new_list in range(len(centres))}
    for data in dataset:
        clusters[distances(data, centres)[1]].append(data)
    return clusters


def new_centres(clusters):
    """
    

    Parameters
    ----------
    clusters : TYPE
        DESCRIPTION.

    Returns
    -------
    centres : TYPE
        DESCRIPTION.

    """
    centres = []
    for i in clusters:
        centres.append(gravity(clusters[i]))
    return centres


def lloyd(dataset, k):
    """
    

    Parameters
    ----------
    dataset : TYPE
        DESCRIPTION.
    k : TYPE
        DESCRIPTION.

    Returns
    -------
    centres : TYPE
        DESCRIPTION.

    """
    centres = dataset[:k]
    count = 0
    while count < 100:
        clusters = centres_to_clusters(dataset, centres)
        centres = new_centres(clusters)
        count += 1
    return centres


# =============================================================================
# SECTION 2
# =============================================================================
def gravity_centres(dataset, hidden_vector):
    import numpy as np
    denom = sum(hidden_vector)
    num = []
    for i, data in enumerate(dataset):
        value = np.dot(hidden_vector[i], data)
        num.append(value)
    return sum(num)/denom


# =============================================================================
# Testing
# =============================================================================
if __name__ == "__main__":
    my.flash()
    
    
    