function [Stiffness, Impedance, VolumeFraction] = Cube(NameData, Grid_Number, cube_size, u_disp, rho, E, nu)
    
    % Directory_files='/home/rkm41/GraphDataZ/';
    Directory_files='C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\';
    
    import com.comsol.model.*
    import com.comsol.model.util.*
    
    Grid_Number = Grid_Number-1;

    name_of_model=[NameData,'.mph'];
    
    model = ModelUtil.create('Model');

    name_of_model=[NameData,'.mph'];
    nodes=load([Directory_files,NameData,'.mat']);
    Number_of_Cubes=length(nodes.cubes);
    CenterNodes=zeros(Number_of_Cubes,3);
    CenterNodes(:,:)=nodes.cubes(:,:);

    model.component.create('comp1', true);
       
    model.param.set('CellSize', [num2str(cube_size), ' [mm]']);
    model.param.set('SmallCubesSize', [num2str(cube_size/Grid_Number), ' [mm]']);
    model.param.set('u_disp', [num2str(u_disp), ' [mm]']);
    
    model.component('comp1').geom.create('geom1', 3);
    model.component('comp1').geom('geom1').lengthUnit('mm');
    % model.component('comp1').geom('geom1').geomRep('comsol');
    model.component('comp1').geom('geom1').geomRep('cadps');
    
    for i=1:Number_of_Cubes
        
        model.component('comp1').geom('geom1').create(['blk',num2str(i)], 'Block');
        model.component('comp1').geom('geom1').feature(['blk',num2str(i)]).set('size', {'SmallCubesSize' 'SmallCubesSize' 'SmallCubesSize'});
        model.component('comp1').geom('geom1').feature(['blk',num2str(i)]).set('base', 'center');
        model.component('comp1').geom('geom1').feature(['blk',num2str(i)]).set('pos', CenterNodes(i,:));
        
    end
    
    model.component('comp1').geom('geom1').run;
    
    model.component('comp1').material.create('mat1', 'Common');
    model.component('comp1').material('mat1').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');

    model.component('comp1').massProp.create('mass1', 'MassProperties');
    model.component('comp1').massProp('mass1').selection.geom('geom1', 3);
    model.component('comp1').massProp('mass1').selection.all;
    model.component('comp1').massProp('mass1').set('densitySource', 'fromPhysics');

    model.component('comp1').material('mat1').propertyGroup('def').set('density', num2str(rho));
    model.component('comp1').material('mat1').propertyGroup('Enu').set('E', [num2str(E), ' [GPa]']);
    model.component('comp1').material('mat1').propertyGroup('Enu').set('nu', num2str(nu));
    
    % model.save([Directory_files+name_of_model])

    model.component('comp1').selection.create('box1', 'Box');
    model.component('comp1').selection('box1').set('entitydim', 2);
    model.component('comp1').selection('box1').label('XP');
    model.component('comp1').selection('box1').set('xmin', '0.9999*CellSize');
    model.component('comp1').selection('box1').set('xmax', '1.0001*CellSize');
    model.component('comp1').selection('box1').set('ymin', 0);
    model.component('comp1').selection('box1').set('ymax', 'CellSize');
    model.component('comp1').selection('box1').set('zmin', 0);
    model.component('comp1').selection('box1').set('zmax', 'CellSize');
    model.component('comp1').selection('box1').set('condition', 'inside');
    
    model.component('comp1').selection.create('box2', 'Box');
    model.component('comp1').selection('box2').set('entitydim', 2);
    model.component('comp1').selection('box2').label('XN');
    model.component('comp1').selection('box2').set('xmin', '-0.0001');
    model.component('comp1').selection('box2').set('xmax', '0.0001');
    model.component('comp1').selection('box2').set('ymin', 0);
    model.component('comp1').selection('box2').set('ymax', 'CellSize');
    model.component('comp1').selection('box2').set('zmin', 0);
    model.component('comp1').selection('box2').set('zmax', 'CellSize');
    model.component('comp1').selection('box2').set('condition', 'inside');
    
    model.component('comp1').selection.create('box3', 'Box');
    model.component('comp1').selection('box3').set('entitydim', 2);
    model.component('comp1').selection('box3').label('YP');
    model.component('comp1').selection('box3').set('xmin', 0);
    model.component('comp1').selection('box3').set('xmax', 'CellSize');
    model.component('comp1').selection('box3').set('ymin', '0.9999*CellSize');
    model.component('comp1').selection('box3').set('ymax', '1.0001*CellSize');
    model.component('comp1').selection('box3').set('zmin', 0);
    model.component('comp1').selection('box3').set('zmax', 'CellSize');
    model.component('comp1').selection('box3').set('condition', 'inside');
    
    model.component('comp1').selection.create('box4', 'Box');
    model.component('comp1').selection('box4').set('entitydim', 2);
    model.component('comp1').selection('box4').label('YN');
    model.component('comp1').selection('box4').set('xmin', 0);
    model.component('comp1').selection('box4').set('xmax', 'CellSize');
    model.component('comp1').selection('box4').set('ymin', '-0.0001');
    model.component('comp1').selection('box4').set('ymax', '0.0001');
    model.component('comp1').selection('box4').set('zmin', 0);
    model.component('comp1').selection('box4').set('zmax', 'CellSize');
    model.component('comp1').selection('box4').set('condition', 'inside');
    
    model.component('comp1').selection.create('box5', 'Box');
    model.component('comp1').selection('box5').set('entitydim', 2);
    model.component('comp1').selection('box5').label('ZP');
    model.component('comp1').selection('box5').set('xmin', 0);
    model.component('comp1').selection('box5').set('xmax', 'CellSize');
    model.component('comp1').selection('box5').set('ymin', 0);
    model.component('comp1').selection('box5').set('ymax', 'CellSize');
    model.component('comp1').selection('box5').set('zmin', '0.9999*CellSize');
    model.component('comp1').selection('box5').set('zmax', '1.0001*CellSize');
    model.component('comp1').selection('box5').set('condition', 'inside');
    
    model.component('comp1').selection.create('box6', 'Box');
    model.component('comp1').selection('box6').set('entitydim', 2);
    model.component('comp1').selection('box6').label('ZN');
    model.component('comp1').selection('box6').set('xmin', 0);
    model.component('comp1').selection('box6').set('xmax', 'CellSize');
    model.component('comp1').selection('box6').set('ymin', 0);
    model.component('comp1').selection('box6').set('ymax', 'CellSize');
    model.component('comp1').selection('box6').set('zmin', '-0.0001');
    model.component('comp1').selection('box6').set('zmax', '0.0001');
    model.component('comp1').selection('box6').set('condition', 'inside');  
    
    % model.save(model_filename)
    
    model.component('comp1').physics.create('solid', 'SolidMechanics', 'geom1');

    model.component('comp1').physics('solid').prop('ShapeProperty').set('order_displacement', 1);

    model.component('comp1').physics('solid').create('disp1', 'Displacement2', 2);
    model.component('comp1').physics('solid').feature('disp1').selection.named('box1');
    model.component('comp1').physics('solid').feature('disp1').set('Direction', [1; 1; 1]);
    model.component('comp1').physics('solid').feature('disp1').set('U0', {'u_disp'; '0'; '0'});
    
    model.component('comp1').physics('solid').create('disp2', 'Displacement2', 2);
    model.component('comp1').physics('solid').feature('disp2').selection.named('box2');
    model.component('comp1').physics('solid').feature('disp2').set('Direction', [1; 1; 1]);
    model.component('comp1').physics('solid').feature('disp2').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('solid').create('disp3', 'Displacement2', 2);
    model.component('comp1').physics('solid').feature('disp3').selection.named('box3');
    model.component('comp1').physics('solid').feature('disp3').set('Direction', [0; 1; 1]);
    model.component('comp1').physics('solid').feature('disp3').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('solid').create('disp4', 'Displacement2', 2);
    model.component('comp1').physics('solid').feature('disp4').selection.named('box4');
    model.component('comp1').physics('solid').feature('disp4').set('Direction', [0; 1; 1]);
    model.component('comp1').physics('solid').feature('disp4').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('solid').create('disp5', 'Displacement2', 2);
    model.component('comp1').physics('solid').feature('disp5').selection.named('box5');
    model.component('comp1').physics('solid').feature('disp5').set('Direction', [0; 1; 1]);
    model.component('comp1').physics('solid').feature('disp5').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('solid').create('disp6', 'Displacement2', 2);
    model.component('comp1').physics('solid').feature('disp6').selection.named('box6');
    model.component('comp1').physics('solid').feature('disp6').set('Direction', [0; 1; 1]);
    model.component('comp1').physics('solid').feature('disp6').set('U0', {'0'; '0'; '0'});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    model.component('comp1').physics('solid').create('disp7', 'Displacement2', 2);
    model.component('comp1').physics('solid').feature('disp7').selection.named('box1');
    model.component('comp1').physics('solid').feature('disp7').set('Direction', [1; 0; 1]);
    model.component('comp1').physics('solid').feature('disp7').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('solid').create('disp8', 'Displacement2', 2);
    model.component('comp1').physics('solid').feature('disp8').selection.named('box2');
    model.component('comp1').physics('solid').feature('disp8').set('Direction', [1; 0; 1]);
    model.component('comp1').physics('solid').feature('disp8').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('solid').create('disp9', 'Displacement2', 2);
    model.component('comp1').physics('solid').feature('disp9').selection.named('box3');
    model.component('comp1').physics('solid').feature('disp9').set('Direction', [1; 1; 1]);
    model.component('comp1').physics('solid').feature('disp9').set('U0', {'0'; 'u_disp'; '0'});
    
    model.component('comp1').physics('solid').create('disp10', 'Displacement2', 2);
    model.component('comp1').physics('solid').feature('disp10').selection.named('box4');
    model.component('comp1').physics('solid').feature('disp10').set('Direction', [1; 1; 1]);
    model.component('comp1').physics('solid').feature('disp10').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('solid').create('disp11', 'Displacement2', 2);
    model.component('comp1').physics('solid').feature('disp11').selection.named('box5');
    model.component('comp1').physics('solid').feature('disp11').set('Direction', [1; 0; 1]);
    model.component('comp1').physics('solid').feature('disp11').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('solid').create('disp12', 'Displacement2', 2);
    model.component('comp1').physics('solid').feature('disp12').selection.named('box6');
    model.component('comp1').physics('solid').feature('disp12').set('Direction', [1; 0; 1]);
    model.component('comp1').physics('solid').feature('disp12').set('U0', {'0'; '0'; '0'});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    model.component('comp1').physics('solid').create('disp13', 'Displacement2', 2);
    model.component('comp1').physics('solid').feature('disp13').selection.named('box1');
    model.component('comp1').physics('solid').feature('disp13').set('Direction', [1; 1; 0]);
    model.component('comp1').physics('solid').feature('disp13').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('solid').create('disp14', 'Displacement2', 2);
    model.component('comp1').physics('solid').feature('disp14').selection.named('box2');
    model.component('comp1').physics('solid').feature('disp14').set('Direction', [1; 1; 0]);
    model.component('comp1').physics('solid').feature('disp14').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('solid').create('disp15', 'Displacement2', 2);
    model.component('comp1').physics('solid').feature('disp15').selection.named('box3');
    model.component('comp1').physics('solid').feature('disp15').set('Direction', [1; 1; 0]);
    model.component('comp1').physics('solid').feature('disp15').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('solid').create('disp16', 'Displacement2', 2);
    model.component('comp1').physics('solid').feature('disp16').selection.named('box4');
    model.component('comp1').physics('solid').feature('disp16').set('Direction', [1; 1; 0]);
    model.component('comp1').physics('solid').feature('disp16').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('solid').create('disp17', 'Displacement2', 2);
    model.component('comp1').physics('solid').feature('disp17').selection.named('box5');
    model.component('comp1').physics('solid').feature('disp17').set('Direction', [1; 1; 1]);
    model.component('comp1').physics('solid').feature('disp17').set('U0', {'0'; '0'; 'u_disp'});
    
    model.component('comp1').physics('solid').create('disp18', 'Displacement2', 2);
    model.component('comp1').physics('solid').feature('disp18').selection.named('box6');
    model.component('comp1').physics('solid').feature('disp18').set('Direction', [1; 1; 1]);
    model.component('comp1').physics('solid').feature('disp18').set('U0', {'0'; '0'; '0'});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % model.save([Directory_files,name_of_model])

    model.component('comp1').mesh.create('mesh1');
    model.component('comp1').mesh('mesh1').autoMeshSize(4);
    model.component('comp1').mesh('mesh1').run;

    % model.save([Directory_files,name_of_model])

    model.study.create('std1');
    model.study('std1').create('stat', 'Stationary');
    model.study('std1').feature('stat').set('useadvanceddisable', true);
    model.study('std1').feature('stat').set('disabledphysics', {'solid/disp7' 'solid/disp8' 'solid/disp9' 'solid/disp10' 'solid/disp11' 'solid/disp12' 'solid/disp13' 'solid/disp14' 'solid/disp15' 'solid/disp16'  ...
    'solid/disp17' 'solid/disp18'});
    
    model.study.create('std2');
    model.study('std2').create('stat', 'Stationary');
    model.study('std2').feature('stat').set('useadvanceddisable', true);
    model.study('std2').feature('stat').set('disabledphysics', {'solid/disp1' 'solid/disp2' 'solid/disp3' 'solid/disp4' 'solid/disp5' 'solid/disp6' 'solid/disp13' 'solid/disp14' 'solid/disp15' 'solid/disp16'  ...
    'solid/disp17' 'solid/disp18'});
    
    model.study.create('std3');
    model.study('std3').create('stat', 'Stationary');
    model.study('std3').feature('stat').set('useadvanceddisable', true);
    model.study('std3').feature('stat').set('disabledphysics', {'solid/disp1' 'solid/disp2' 'solid/disp3' 'solid/disp4' 'solid/disp5' 'solid/disp6' 'solid/disp7' 'solid/disp8' 'solid/disp9' 'solid/disp10'  ...
    'solid/disp11' 'solid/disp12'});
    
    model.save([Directory_files,name_of_model])
    
    model.study('std1').run;
    model.study('std2').run;
    model.study('std3').run;
    
    % model.save([Directory_files,name_of_model])
    
    UEnergy1=mphglobal(model,'solid.Ws_tot', 'dataset', 'dset1');
    UEnergy2=mphglobal(model,'solid.Ws_tot', 'dataset', 'dset2');
    UEnergy3=mphglobal(model,'solid.Ws_tot', 'dataset', 'dset3');
    
    C11=2*UEnergy1/(cube_size/1000)^3/(u_disp/cube_size)^2;
    C22=2*UEnergy2/(cube_size/1000)^3/(u_disp/cube_size)^2;
    C33=2*UEnergy3/(cube_size/1000)^3/(u_disp/cube_size)^2;
    
    MassValue=mphglobal(model,'mass1.mass');
    DensityStructure=MassValue/(cube_size/1000)^3;
    
    Impedance1 = sqrt(C11*DensityStructure);
    Impedance2 = sqrt(C22*DensityStructure);
    Impedance3 = sqrt(C33*DensityStructure);
    
    VolumeFraction = MassValue/(rho*(cube_size/1000)^3);
    
    Stiffness = [C11, C22, C33];
    Impedance = [Impedance1, Impedance2, Impedance3];
    
    output = struct();
    output.Impedance = Impedance;
    output.Stiffness = Stiffness;
    output.VolumeFraction = VolumeFraction;
    output.DensityStructure = DensityStructure;
    
    % Save the new parameters to the existing .mat file, appending the data
    save([Directory_files, NameData, '.mat'], '-struct', 'output', '-append');
    
end