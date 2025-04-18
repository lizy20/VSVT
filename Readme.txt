
% % % MIT License
% % % Copyright (c) 2025 [Zeyang Li,M2O group, School of engineering, Cardiff University]

% In order to reduce code redundancy,  some functionality has been removed or annotated with interface.  
% Interested readers can debug the codes by themselves according to the requirement，
         refer to the paperwork{https://doi.org/10.1016/j.cma.2024.117378}.
% For the sake of computation efficiency on normal desktop, we certainly 
        adjust serval parameters in the codes, with decreasing  volume//compliance accuracy in 
        de-homogenization process.  Please refer to the comment in codes for reference.

Any other questions about the codes， please contact author:   www.m2olab.uk// 4692420536@qq.com

% 【【【 For citation of this code,  please refer to https://doi.org/10.1016/j.cma.2024.117378 】】】

%-----------------------------------------------------------------------------------------------------------------------
======【2D Programmable, anisotropic bio-inspired porous material optimization method by varied-shape Voronoi tessellation】===================

A MATLAB App Designer-based graphical interface for Structural optimization.
## Requirements
-  MATLAB R2023b or newer
-  attached open-access  ‘Poisson disc sampling’， ‘SlanCM colormap’.

-  

## 🚀 Getting Started 

   1.run 'run_main_app.mlapp'('run_main_app.m') file for an interactive GUI with App Designer.
   
   2.Given the info of the typical rectangle design domain and constraint

	- Rectangle Domain size, meshgrid number, global volume fraction constraint.

	- Boundary condition { cantilever， MBB， Unique-direction compression} refer to the paper.

	  Typically, for Uni-direction compression case, we used extra element direction constraint, alpha = 89~90°，to achieve
                 results in paper {https://doi.org/10.1016/j.cma.2024.117378}, marked as case-4 compress(paper)
	   A normal optimization considering {alpha =-90 ~ 90}, please comments out the codes {Row 152~154 in /multito/continous}. However, the performance
	    of the final result may suggest certain decrease, due to the influence by local optimal during the MMA optimization process.


	- Model overview can be check in the Figure 10. Each green nodal force marker correspond 
					normalized force of {1 value, -Y axis direction}.

	- non-structure mesh config can be calculated by User-defined interface. Two cases {demo-1,demo-2} are provided for the reference.

   3. click 'Start' for optimization convergence
	- the convergence history and optimized design variables can be monitored refer to the UI windows. 

   4. After optimization completion, click 'generate model' for building detailed porous geometry.
   
   5. Afterward, click 'generate STL file' for building detailed porous geometry. (The step maybe high time-costing.)

   6. For simple problem, the 'postfunction' processing in 'continous.m' may be non-necessary and can be commented out.
      //Otherwise, adjust R_{min} parameter in {Row26,'/multito/postfunction.m' } for eliminating possible blank interface issue for multi-material optimization.
        R_{min} ≈ 1.15*mesh element length, is usually reasonable.
	
  7.  For sake of converience, we provide an {Aniso-iso switch} to respectively achieve anisotropic//isotropic Voronoi porous structure optmization. 


       
## 🚀 Advanced using
 A series of tuning interface is provided for experienced user in practical application.
    For read-friendly， the basic mesh assignment framework is refer to "An efficient 3D topology optimization code written in Matlab" (Liu K,2013)


    1.Ones can alter and then use 'run_main.m' as main codes without GUI interface.
       - run_main.m     main code for prepared info
       - continuous.m   main optmizaiton code
       - postfunction.m   post-processing code for STL building (non-necessary for simple question)
       
       - run_geo_main.m    dehomogenized the optimization result of {run_main.m}
       - run_build_stl.m     generate stand stl file for designed porous structure.    
    
    2. In GUI interface, ones can use 'user-defined' selection in 'Boudary Condition' list for own mesh-info design.

		---  put a FEMmeshdata.mat file containing {nodes information, element information, 
		     boundary fixed freedom, loading nodes freedom ...}  under '...\meshdata\' directory.

		--- 'continous.m', Row70~93 provides specific processing step for the FEMmeshdata.m. 

			 Any modification can be used here, to provide necessary info, including
			 opt.nele: element number 
			 opt.Node: node coordinates (nx2 matrix)  
			 opt.Element:  cell vector, i^th cell is of [Node id list] for i^th element {nele x 4| nele x3 cell vector}， 
					Meshed elements should be of similar size for later robust geometry computation.  
			 opt.ndof: total nodal DOF  (2n integer)

			 fixednid: list of fixed node id in model,  like[1,2,5,6], suggesting node 1,2,5,6 have limited displacement BC.
			 fixeddofs: list of fixed nodal DOF in model, 2*i-1, 2i corresponds -x,-y direction DOF of i^th node.
				like [1,2, 3, 9, 12], suggesting {-x and -y direction of node 1; -x of node 2,5	;-y of node 6} are fixed in displacement.
						
		         loadingnid: list of loading node id in model;
			 loadingdofs:  list of loading nodal DOF in model;
			 F value: The external force values with respect to loadingdofs.
		 	【For easier mesh building, a meshdata pre-porcessing code named 'User_defined_meshdata_building.m' is provided, focusing on the tranforming
			 of ABAQUS inp model into '.mat'.】

    3. The geometry dehomogenzation process is divided by the Optimization Computation for the flexbiility, as‘’ Compute geometry” button.
         Code framwork as follows:
	 --‘’workflow1.m‘’ ： Dehomogenzation rho--> desired  seed density,{ N_s} and  truss thickness {t} for TO unit.   
			keep N_s =  20 <-->   120 range, for {high computation efficiency} <--> {high compliance accuracy} of homogenized geometry model
			   N_s =60 recommended for realistic  application

	-- "workflow2.m" Varible-density Poisson disk sample depends on N_s;
	--"workflow3.m"  Numerical Voronoi tessellation computation. 
 	--"workflow4.m"  Convexhull & Voronoi truss info computation. 
                 --“workflow5.m”   Build final geometry polyshape
 	--“workflow6.m”    Boundary simplify for final geometry;
	For the sake of convenience, each workflowXX.m will automatically read the former result geoXX.mat and save the processed result as geoXX+1.m//
	So, after the one complete workflow, ones can only adjust single .m file and individually run it for the desired effect, by using  GUI- {All steps} selection panel.
	This may be useful for testing different boudanry thickness and  simplification. 


   
