

```R
require(Lightbulb)
```

    Loading required package: Lightbulb
    Loading required package: data.table
    Loading required package: matrixStats
    Loading required package: Matrix
    Loading required package: gplots
    
    Attaching package: ‘gplots’
    
    The following object is masked from ‘package:stats’:
    
        lowess
    
    Loading required package: ggplot2
    Loading required package: Seurat
    Loading required package: cowplot
    
    Attaching package: ‘cowplot’
    
    The following object is masked from ‘package:ggplot2’:
    
        ggsave
    
    
    Attaching package: ‘Lightbulb’
    
    The following objects are masked from ‘package:matrixStats’:
    
        colMaxs, colMins, colSds, rowMaxs, rowMins, rowSds
    


## loading data


```R
Exp_Seurat=readRDS("Lightbulb_example_ExpSeurat.rds")
```

### Super cell creation

Exp_Seurat is a Seurat object. Seurat object stores the log2(TPM) matrix in @data slot. Our super cell function can take a seurat object as input or a cell-gene matrix as input. For details about generating your own Seurat object, please refer to Seurat manual (https://satijalab.org/seurat/get_started.html)


```R
Exp_Seurat@data[1:5,1:5]
```


    5 x 5 sparse Matrix of class "dgCMatrix"
            CAGTCCTTCCAAGTAC-1 TAAGTGCTCTCTGAGA-1 CGCTTCATCTTACCTA-1
    Mrpl15            .                  8.470585                  .
    Lypla1            .                  8.470585                  .
    Tcea1             8.633139           8.470585                  .
    Atp6v1h           .                  .                         .
    Rb1cc1            .                  8.470585                  .
            TCTTCGGAGTTTAGGA-1 AATCCAGCATAGACTC-1
    Mrpl15                   .                  .
    Lypla1                   .                  .
    Tcea1                    .                  .
    Atp6v1h                  .                  .
    Rb1cc1                   .                  .



```R
Exp_Seurat@meta.data[1:10,]
```


<table>
<thead><tr><th></th><th scope=col>nGene</th><th scope=col>nUMI</th><th scope=col>orig.ident</th><th scope=col>batch</th><th scope=col>res.1</th><th scope=col>tissue</th><th scope=col>timepoint</th><th scope=col>test</th><th scope=col>res.3</th></tr></thead>
<tbody>
	<tr><th scope=row>CAGTCCTTCCAAGTAC-1</th><td>945          </td><td>2525         </td><td>SeuratProject</td><td>1            </td><td>8            </td><td>Spleen       </td><td>d0           </td><td>0            </td><td>0            </td></tr>
	<tr><th scope=row>TAAGTGCTCTCTGAGA-1</th><td>891          </td><td>2827         </td><td>SeuratProject</td><td>1            </td><td>8            </td><td>Spleen       </td><td>d0           </td><td>0            </td><td>0            </td></tr>
	<tr><th scope=row>CGCTTCATCTTACCTA-1</th><td>981          </td><td>3012         </td><td>SeuratProject</td><td>1            </td><td>8            </td><td>Spleen       </td><td>d0           </td><td>0            </td><td>0            </td></tr>
	<tr><th scope=row>TCTTCGGAGTTTAGGA-1</th><td>744          </td><td>2212         </td><td>SeuratProject</td><td>1            </td><td>8            </td><td>Spleen       </td><td>d0           </td><td>0            </td><td>0            </td></tr>
	<tr><th scope=row>AATCCAGCATAGACTC-1</th><td>890          </td><td>2455         </td><td>SeuratProject</td><td>1            </td><td>8            </td><td>Spleen       </td><td>d0           </td><td>0            </td><td>0            </td></tr>
	<tr><th scope=row>CTTAGGAGTGTTAAGA-1</th><td>680          </td><td>1628         </td><td>SeuratProject</td><td>1            </td><td>8            </td><td>Spleen       </td><td>d0           </td><td>0            </td><td>0            </td></tr>
	<tr><th scope=row>TGACTAGTCCTTTCGG-1</th><td>996          </td><td>3462         </td><td>SeuratProject</td><td>1            </td><td>8            </td><td>Spleen       </td><td>d0           </td><td>0            </td><td>0            </td></tr>
	<tr><th scope=row>GACACGCGTACCGTAT-1</th><td>717          </td><td>2552         </td><td>SeuratProject</td><td>1            </td><td>8            </td><td>Spleen       </td><td>d0           </td><td>0            </td><td>0            </td></tr>
	<tr><th scope=row>CTAATGGCACTGTGTA-1</th><td>716          </td><td>1631         </td><td>SeuratProject</td><td>1            </td><td>8            </td><td>Spleen       </td><td>d0           </td><td>0            </td><td>0            </td></tr>
	<tr><th scope=row>TTAGGCATCCTGCCAT-1</th><td>856          </td><td>2483         </td><td>SeuratProject</td><td>1            </td><td>8            </td><td>Spleen       </td><td>d0           </td><td>0            </td><td>0            </td></tr>
</tbody>
</table>



#### using seurat object as input
When taking Seurat object as input, the calculation is based on @scale.data slot. We need to make sure that @scale.data is addible. In this example, we simply put the TPM matrix in the @scale.data slot.

*k_merge* an integer indicates the number of cells we want to merge to create a super-cell.  
*n* an integer of the number of super-cell we want to get.  
*sampling_ref* a vector of group_id for each cell. The sampling process will try to make sure the group_id distribution remains the same after sampling.  


```R
Exp_Seurat@scale.data=2^Exp_Seurat@data-1

#calculate super cell
supercell_mat=Super_cell_creation(Exp_Seurat,k_merge = 50,n=6000,sampling_ref = Exp_Seurat@meta.data$batch)
supercell_mat=log2(supercell_mat+1)
```

    Seurat object detected as input
    MNN done: time consumed: 0 hr 1 min 51.37 s
    merging finished: time consumed: 0 hr 4 min 54.86 s


#### using cell gene matrix  as input
The cell gene matrix is a matrix which each row represents a cell and each column is a gene. We need to make sure the value of cell-gene matrix is addible. In this example, we simply use the TPM matrix.


```R
TPM_matrix=as.matrix(t(2^Exp_Seurat@data-1))

#calculate super cell
supercell_mat=Super_cell_creation(TPM_matrix,k_merge = 50,n=6000,sampling_ref = Exp_Seurat@meta.data$batch)
supercell_mat=log2(supercell_mat+1)
```

    Assuming input is cell-gene matrix
    MNN done: time consumed: 0 hr 1 min 34.67 s
    merging finished: time consumed: 0 hr 5 min 26.77 s


The row name of the super-cell matrix indicate the position of the cell center in the original data. For example in the following super-cell matrix, the first super-cell is centered on cell 193 in orginal data (row 193 in TPM_matrix or column 193 in Exp_Seurat@scale.data).


```R
supercell_mat[1:5,1:5]
```


<table>
<thead><tr><th></th><th scope=col>Mrpl15</th><th scope=col>Lypla1</th><th scope=col>Tcea1</th><th scope=col>Atp6v1h</th><th scope=col>Rb1cc1</th></tr></thead>
<tbody>
	<tr><th scope=row>193</th><td>5.458230</td><td>5.265078</td><td>5.926450</td><td>4.828982</td><td>5.057072</td></tr>
	<tr><th scope=row>276</th><td>5.754755</td><td>5.342182</td><td>4.219048</td><td>5.058700</td><td>4.709162</td></tr>
	<tr><th scope=row>1234</th><td>6.467162</td><td>5.142626</td><td>5.424017</td><td>5.212062</td><td>4.988101</td></tr>
	<tr><th scope=row>651</th><td>6.353871</td><td>4.017780</td><td>5.602713</td><td>4.489083</td><td>3.023346</td></tr>
	<tr><th scope=row>495</th><td>5.606721</td><td>5.507785</td><td>5.874736</td><td>4.482352</td><td>3.114476</td></tr>
</tbody>
</table>




```R

```
