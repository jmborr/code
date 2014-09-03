import numpy

def rmsd(crds1, crds2):
      
      """Returns RMSD between 2 sets of [nx3] numpy array"""
      assert(crds1.shape[1]==3)
      assert(crds1.shape==crds2.shape)
      n_vec = numpy.shape(crds1)[0]
      correlation_matrix = numpy.dot(numpy.transpose(crds1), crds2)
      v, s, w_tr = numpy.linalg.svd(correlation_matrix)
      is_reflection = (numpy.linalg.det(v) * numpy.linalg.det(w_tr)) < 0.0
      if is_reflection: s[-1] = - s[-1]
      E0 = sum(sum(crds1 * crds1)) + sum(sum(crds2 * crds2))
      rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
      rmsd_sq = max([rmsd_sq, 0.0])
      return numpy.sqrt(rmsd_sq)

def optimal_superposition(crds1, crds2):
    """Returns best-fit rotation matrix as [3x3] numpy matrix"""
    assert(crds1.shape[1] == 3)
    assert(crds1.shape == crds2.shape)
    correlation_matrix = numpy.dot(numpy.transpose(crds1), crds2)
    v, s, w_tr = numpy.linalg.svd(correlation_matrix)
    is_reflection = (numpy.linalg.det(v) * numpy.linalg.det(w_tr)) < 0.0
    if is_reflection: v[:,-1] = -v[:,-1]
    return numpy.dot(v, w_tr)

                                                
