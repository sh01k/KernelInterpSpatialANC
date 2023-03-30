# KernelInterpSpatialANC

## Description
An active noise control (ANC) method to reduce noise over a region in space based on kernel interpolation of sound field. This repository provides codes for reproducing results in the frequency domain described in the following article.

- S. Koyama, J. Brunnström, H. Ito, N. Ueno, and H. Saruwatari, "Spatial Active Noise Control Based on Kernel Interpolation of Sound Field," *IEEE/ACM Transactions on Audio, Speech, and Language Processing*, DOI: 10.1109/TASLP.2021.3107983, 2021.

The article is open access on [IEEE Xplore](https://doi.org/10.1109/TASLP.2021.3107983).

### Abstract
An active noise control (ANC) method to reduce noise over a region in space based on kernel interpolation of sound field is proposed. Current methods of spatial ANC are largely based on spherical or circular harmonic expansion of the sound field, where the geometry of the error microphone array is restricted to a simple one such as a sphere or circle. We instead apply the kernel interpolation method, which allows for the estimation of a sound field in a continuous region with flexible array configurations. The interpolation scheme is used to derive adaptive filtering algorithms for minimizing the acoustic potential energy inside a target region. A practical time-domain algorithm is also developed together with its computationally efficient block-based equivalent. We conduct experiments to investigate the achievable level of noise reduction in a two-dimensional free space, as well as adaptive broadband noise control in a three-dimensional reverberant space.  The experimental results indicated that the proposed method outperforms the multipoint-pressure-control-based method in terms of regional noise reduction. 

## License
[MIT](https://github.com/sh01k/KernelInterpSpatialANC/blob/main/LICENSE)

## Author
- [Shoichi Koyama](https://www.sh01.org) 
- Jesper Brunnström
- Hayato Ito
- [Natsuki Ueno](https://natsuenono.github.io/)
- [Hiroshi Saruwatari](https://researchmap.jp/read0102891/)