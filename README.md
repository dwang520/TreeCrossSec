# TreeCrossSec
Tree stem cross-section reconstruction with free-form curves from point clouds
# Usage
VV = TreeCrossSec(Sp, res_thres, Forder) <br/>

Input:<br/>
Sp: nxm point cloud (unit in m). The dimension m should be at least 2.<br/>
optional:<br/>
res_thres: residual tolerance e.g., 0.01 (cm)<br/>
Forder: Fourier series order (1 - 8)<br/>

Output:<br/>
VV: a structure variable contains diameter estimations from various methods, and reconstructed points in ModelX&ModelY. Usually FouCirc gives the best diameter estimation.

# Examples
Example 1 <br/>
VV = TreeCrossSec(Sp, 0.01, 8);
![example 1](e1.png)
Example 2 <br/> 
(this method only works perfectly for convex shapes) <br/>
VV = TreeCrossSec(Sp, 0.05, 8);
![example 2](e2.png)
Example 3 <br/>
VV = TreeCrossSec(Sp, 0.01, 6);
![example 3](e3.png)
Example 4 <br/>
VV = TreeCrossSec(Sp, 0.01, 6);
![example 4](e4.png)
# Bibtex
@article{xxxx, <br/>
  title={Reconstructing stem cross section shapes from terrestrial laser scanning}, <br/>
  author={Wang, Di and Kankare, Ville and Puttonen, Eetu and Hollaus, Markus and Pfeifer, Norbert}, <br/>
  journal={IEEE Geoscience and Remote Sensing Letters}, <br/>
  volume={14}, <br/>
  number={2}, <br/>
  pages={272--276}, <br/>
  year={2017}, <br/>
  publisher={IEEE} <br/>
}
# Contact
Di Wang <br/>
di.wang@aalto.fi
