# NPS-analysis-JOVE

### mechano-NPS data analysis for the JOVE protocol paper
*copyright 2022 Kristen Cotner, Brian Li, Alan Dong*

Kim, J., Han, S., Lei, A. *et al.* Characterizing cellular mechanical phenotypes with mechano-node-pore sensing. *Microsyst Nanoeng* **4**, 17091 (2018). https://doi.org/10.1038/micronano.2017.91

#### example: data processing
`>> output_table = mNPS_procJOVE(filepath, ch_height, De_np, wC, thresholds, sampleRate);`  
`>> disp(output_table);`  
`>> summary(output_table);`  
`>> histogram(output_table.wCDI);`  
`>> output_array = output_table.cell_data;`  
`>> writetable(output_table, 'output_file.txt');`

#### example: quick data visualization
`>> [y_smoothed, y_downsampled, ~] = mNPS_fastQC(data)`  
`>> plot(y_downsampled)`

### output of `mNPS_procJOVE()`
* `output_table` = _n_ x 19 table, where _n_ is the number of processed cell events
	* units of each column: `output_table.Properties.VariableUnits`
* `output_table.diameter` = cell diameter [µm]
* `output_table.def_diameter` = transverse diameter of the cell in the contraction segment [µm]
	* transverse deformation is defined as `def_diameter/diameter`
* `output_table.wCDI` = whole-cell deformability index
	* dimensionless measure of resistance to deformation
	* _wCDI_ has an inverse relationship to Young's modulus
* `output_table.rec_cat` = the first segment where the cell was fully recovered after deformation [categorical]
	* "recovery" is defined as the current pulse returning to its pre-deformation amplitude, within 8% tolerance
	* `0` = the cell was already recovered by the first recovery segment
	* `1` = the cell did not recover until the second recovery segment
	* `2` = the cell did not recover until the third recovery segment
	* `Inf` = the cell did not recover within the NPS channel
* `output_table.rec_time` = time for the cell to recover after deformation [ms]
	* elapsed time between the cell exiting the contraction segment and entering the first segment where it was recovered
	* "recovery" is defined as the current pulse returning to its pre-deformation amplitude, within 8% tolerance
	* although this is reported as a continuous value, it is discretized for each cell by the three recovery segments
	* cells that did not recover within the NPS channel have a recovery time of `Inf`

#### example data:

* 20220617_A549_dev2B_w12_p25_try1.mat
	* acquired by RR on 06.17.2022 using the main SRL NPS platform
	* devices: sNPS_ver2.1 (wC=12), channel height ~= 30µm
	* A549 cells, passage 3, in PBS + 2% FBS

* 20220616_BEAS2B_Dev1B_w10_p30_try1.mat
	* acquired by AL on 06.16.2022 using the main SRL NPS platform
	* devices: sNPS_ver2.1 (wC=10), channel height ~= 30µm
	* BEAS2B cells

* 20220616_BEAS2B_Dev1B_w10_p30_try1_crop.mat
	* second half of data from 20220616_BEAS2B_Dev1B_w10_p30_try1.mat

#### example results:

* output_j12.mat
	* results of running `mNPS_procJOVE()` using data from 20220617_A549_dev2B_w12_p25_try1.mat (using testscript.m)
	* includes source data file names; processing parameters & default thresholds; and output table

* output_j10c.mat
	* results of running `mNPS_procJOVE()` using data from 20220616_BEAS2B_Dev1B_w10_p30_try1_crop.mat (using testscript.m)
	* includes source data file names; processing parameters & default thresholds; and output table
