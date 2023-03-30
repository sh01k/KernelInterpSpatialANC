import numpy as np
import scipy.spatial.distance as distfuncs
import scipy.special as special

def uniformFilledRectangle(numPoints, lim):
    if len(lim) == 2:
        x = np.linspace(lim[0], lim[1], numPoints[0])
        y = np.linspace(lim[0], lim[1], numPoints[1])
    elif len(lim) == 4:
        x = np.linspace(lim[0], lim[2], numPoints[0])
        y = np.linspace(lim[1], lim[3], numPoints[1])
    
    [xGrid, yGrid] = np.meshgrid(x, y)
    points = np.vstack((xGrid.flatten(), yGrid.flatten())).T

    return points

def equidistantRectangle(numPoints, dims, offset=0.5):
    if numPoints == 0:
        return np.zeros((0, 2))
    totalLength = 2 * (dims[0] + dims[1])
    pointDist = totalLength / numPoints

    points = np.zeros((numPoints, 2))
    if numPoints < 4:
        points = equidistantRectangle(4, dims)
        pointChoices = np.random.choice(4, numPoints, replace=False)
        points = points[pointChoices, :]
    else:
        lengths = [dims[0], dims[1], dims[0], dims[1]]
        xVal = [-dims[0] / 2, dims[0] / 2, dims[0] / 2, -dims[0] / 2]
        yVal = [-dims[1] / 2, -dims[1] / 2, dims[1] / 2, dims[1] / 2]

        startPos = pointDist * offset
        xFac = [1, 0, -1, 0]
        yFac = [0, 1, 0, -1]
        numCounter = 0

        for i in range(4):
            numAxisPoints = 1 + int((lengths[i] - startPos) / pointDist)
            axisPoints = startPos + np.arange(numAxisPoints) * pointDist
            distLeft = lengths[i] - axisPoints[-1]
            points[numCounter:numCounter + numAxisPoints, 0] = xVal[i] + xFac[i] * axisPoints
            points[numCounter:numCounter + numAxisPoints, 1] = yVal[i] + yFac[i] * axisPoints
            numCounter += numAxisPoints
            startPos = pointDist - distLeft

    return points

def equiangularCircle(numPoints, rad, offset=0):
    ang = np.linspace(0, 2*np.pi, numPoints) + offset*(2*np.pi/numPoints)
    points = np.vstack([rad*np.cos(ang), rad*np.sin(ang)]).T
    return points

def PointSource(k, posSrc, posMic):
    r = distfuncs.cdist(posMic, posSrc)[None,:,:]
    p = -(1j/4) * special.hankel2(0, k[:,None,None] * r)
    return p

###### Kernel interpolation functions ######
def block_rect(rectDims, rng=None):
    if rng is None:
        rng = np.random.RandomState()
    totVol = rectDims[0] * rectDims[1]

    def pointGenerator(numSamples):
        x = rng.uniform(-rectDims[0] / 2, rectDims[0] / 2, numSamples)
        y = rng.uniform(-rectDims[1] / 2, rectDims[1] / 2, numSamples)
        points = np.stack((x, y))
        return points.T

    return pointGenerator, totVol

def block_circ(rad, rng=None):
    if rng is None:
        rng = np.random.RandomState()
    totVol = np.pi * rad**2

    def pointGenerator(numSamples):
        r = np.sqrt(2 * rng.uniform(0, 0.5*rad**2, numSamples))
        ang = rng.uniform(0, 2*np.pi, numSamples)
        points = np.stack((r*np.cos(ang), r*np.sin(ang)))
        return points.T

    return pointGenerator, totVol

def integrableAFunc(k, posErr):
    def intFunc(r):
        distance = np.transpose(distfuncs.cdist(r, posErr), (1, 0))[None,:,:]
        kappa = special.j0(k[:,None,None] * distance)
        funcVal = kappa[:, :, None, :].conj() * kappa[:, None, :, :]
        return funcVal

    return intFunc

def integrableAwFunc(k, posErr, beta=0, ang=0):
    def intFunc(r):
        r_diff = (np.tile(r[None,:,:], (posErr.shape[0],1,1)) - np.tile(posErr[:,None,:], (1,r.shape[0],1)))[None,:,:,:]
        distance = 1j*np.sqrt((beta*np.cos(ang) + 1j*k[:,None,None]*r_diff[:,:,:,0])**2 + (beta*np.sin(ang) + 1j*k[:,None,None]*r_diff[:,:,:,1])**2)
        #distance = 1j*np.sqrt((beta*np.cos(ang) - 1j*k[:,None,None]*r_diff[:,:,:,0])**2 + (beta*np.sin(ang) - 1j*k[:,None,None]*r_diff[:,:,:,1])**2)
        kappa = special.jn(0, distance)
        funcVal = kappa[:, :, None, :].conj() * kappa[:, None, :, :]
        return funcVal

    return intFunc

def integrate(func, pointGenerator, totNumSamples, totalVolume, numPerIter=50):
    """pointGenerator should return np array, [numSpatialDimensions, numPoints]
    func should return np array [funcDims, numPoints],
    where funcDims can be any number of dimensions (in a multidimensional array sense)"""
    samplesPerIter = numPerIter
    numBlocks = int(np.ceil(totNumSamples / samplesPerIter))
    outDims = np.squeeze(func(pointGenerator(1)), axis=-1).shape
    integralVal = np.zeros(outDims)
    print(
        "Starting monte carlo integration \n",
        "Samples per block: ",
        numPerIter,
        "\nTotal samples: ",
        numBlocks * numPerIter,
    )

    for i in range(numBlocks):
        points = pointGenerator(samplesPerIter)
        fVals = func(points)

        newIntVal = (integralVal * i + np.mean(fVals, axis=-1)) / (i + 1)
        print("Block ", i)

        integralVal = newIntVal
    integralVal *= totalVolume
    print("Finished!!")
    return integralVal


