# Image_processing toolbox (not completed)

This repo provides the implementation of several image processing techniques implemented via matlab.

Relative Paper list:
<table>
  <tr>
    <td>Technique name</td>
    <td>Script(s)</td>    
    <td>Paper id</td>
  </tr>
  <tr>
    <td>Bone shape reconstruction</td>
    <td><code>image_processing/demo/bone_shape_reconstruct_from_dcm.m</code></td>    
    <td><code>Leng2025.05.02.651799</code></td>
  </tr>
</table>

# Demos

# Usage
## Setup
```
[$DOWNLOAD_DIR]/image_processing/           
├── data/[$task_name]/
    
├── image_proc-toolbox/
    ├── ind2lab.m
    ├── keep_remove_ranked_vol.m
    ├── lab2ind.m
    ├── serial_dilate_erode_vol.m
    ├── sub2ind_nd.m
    ├── vol_rot_gif.m
    
├── file_proc-toolbox/
    ├── get_dirs.m
    ├── get_predirs.m
    ├── imfiles2giffiles.m
    
```

# For the citation
The usage of the bone reconstruction script should cite:
```bibtex
@article {Leng2025.05.02.651799,
	author = {Leng, Houfu and jiang, jiahao and Gassner, Katja and Midha, Swati and Justo-Mendez, Raquel and Zheng, Jianqing and Hall, Timothy and Luo, Lin and West, Suzanne D and Vincent, Tonia L. and Wann, Angus and Patel, Kashyap A and Poulton, Joanna and O{\textquoteright}Callaghan, Chris A. and Lechuga-Vieco, Ana Victoria and Simon, Anna Katharina},
	title = {Mitochondrial heterogeneity disrupts osteoclast differentiation and bone resorption by impairing respiratory complex I},
	elocation-id = {2025.05.02.651799},
	year = {2025},
	doi = {10.1101/2025.05.02.651799},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2025/05/02/2025.05.02.651799},
	eprint = {https://www.biorxiv.org/content/early/2025/05/02/2025.05.02.651799.full.pdf},
	journal = {bioRxiv}
}
```
