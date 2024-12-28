function [Stiffness, Impedance, VolumeFraction] = Triangle(NameData, Grid_Number, cube_size, u_disp, rho, E, nu)

    Directory_files='/home/rkm41/GraphDataZ/';
    % Directory_files='C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\';
    
    import com.comsol.model.*
    import com.comsol.model.util.*
    
    name_of_model=[NameData,'.mph'];
    
    model = ModelUtil.create('Model');
    
    thickness=0.05; %mm

    name_of_model=[NameData,'.mph'];
    nodes=load([Directory_files,NameData,'.mat']);

    Number_of_Triangles=length(nodes.triangles);
    CornerNodes=zeros(Number_of_Triangles,9);
    CornerNodes(:,:)=nodes.triangles(:,:);


    model.component.create('comp1', true);
    model.param.set('CellSize', [num2str(cube_size), ' [mm]']);
    model.param.set('u_disp', [num2str(u_disp), ' [mm]']);


    model.component('comp1').geom.create('geom1', 3);
    model.component('comp1').geom('geom1').lengthUnit('mm');
    % model.component('comp1').geom('geom1').geomRep('comsol');
    model.component('comp1').geom('geom1').geomRep('cadps');
    
    for i=1:Number_of_Triangles
        model.component('comp1').geom('geom1').create(['wp',num2str(i)], 'WorkPlane');
        model.component('comp1').geom('geom1').feature(['wp',num2str(i)]).set('planetype', 'coordinates');
        model.component('comp1').geom('geom1').feature(['wp',num2str(i)]).set('genpoints', [CornerNodes(i,1),CornerNodes(i,2),CornerNodes(i,3); CornerNodes(i,4),CornerNodes(i,5),CornerNodes(i,6); CornerNodes(i,7),CornerNodes(i,8),CornerNodes(i,9)]);
        model.component('comp1').geom('geom1').feature(['wp',num2str(i)]).set('unite', true);
        model.component('comp1').geom('geom1').feature(['wp',num2str(i)]).geom.create('pol1', 'Polygon');
        model.component('comp1').geom('geom1').feature(['wp',num2str(i)]).geom.feature('pol1').set('source', 'table');
        Local_Origin=[CornerNodes(i,1);CornerNodes(i,2);CornerNodes(i,3)];
        Local_X=[CornerNodes(i,4);CornerNodes(i,5);CornerNodes(i,6)]-[CornerNodes(i,1);CornerNodes(i,2);CornerNodes(i,3)];
        Local_Y=[CornerNodes(i,7);CornerNodes(i,8);CornerNodes(i,9)]-[CornerNodes(i,1);CornerNodes(i,2);CornerNodes(i,3)];
        Local_Z=cross(Local_X,Local_Y);
        CornerNodes_Global=[CornerNodes(i,1) CornerNodes(i,4) CornerNodes(i,7);CornerNodes(i,2) CornerNodes(i,5) CornerNodes(i,8);CornerNodes(i,3) CornerNodes(i,6) CornerNodes(i,9)];
        CornerNodes_Local=global2localcoord(CornerNodes_Global,'rr',Local_Origin,[Local_X, Local_Y, Local_Z]);
        model.component('comp1').geom('geom1').feature(['wp',num2str(i)]).geom.feature('pol1').set('table', [CornerNodes_Local(1,1), CornerNodes_Local(2,1); CornerNodes_Local(1,2), CornerNodes_Local(2,2);CornerNodes_Local(1,3), CornerNodes_Local(2,3)]);
    end
      
    model.component('comp1').geom('geom1').run;
    model.component('comp1').geom('geom1').run('fin');
    
    model.component('comp1').material.create('mat1', 'Common');
    model.component('comp1').material('mat1').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
    
    model.component('comp1').physics.create('shell', 'Shell', 'geom1');
    
    model.component('comp1').massProp.create('mass1', 'MassProperties');
    model.component('comp1').massProp('mass1').selection.geom('geom1', 2);
    model.component('comp1').massProp('mass1').selection.all;
    model.component('comp1').massProp('mass1').set('densitySource', 'fromPhysics');
    
    model.component('comp1').material('mat1').propertyGroup('def').set('density', num2str(rho));
    model.component('comp1').material('mat1').propertyGroup('Enu').set('E', [num2str(E), ' [GPa]']);
    model.component('comp1').material('mat1').propertyGroup('Enu').set('nu', num2str(nu));

    % model.save([Directory_files,name_of_model])
    
    model.component('comp1').selection.create('box1', 'Box');
    model.component('comp1').selection('box1').set('entitydim', 2);
    model.component('comp1').selection('box1').label('XP');
    model.component('comp1').selection('box1').set('xmin', 'CellSize');
    model.component('comp1').selection('box1').set('xmax', 'CellSize');
    model.component('comp1').selection('box1').set('ymin', 0);
    model.component('comp1').selection('box1').set('ymax', 'CellSize');
    model.component('comp1').selection('box1').set('zmin', 0);
    model.component('comp1').selection('box1').set('zmax', 'CellSize');
    model.component('comp1').selection('box1').set('condition', 'inside');
    
    model.component('comp1').selection.create('box2', 'Box');
    model.component('comp1').selection('box2').set('entitydim', 2);
    model.component('comp1').selection('box2').label('XN');
    model.component('comp1').selection('box2').set('xmin', 0);
    model.component('comp1').selection('box2').set('xmax', 0);
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
    model.component('comp1').selection('box3').set('ymin', 'CellSize');
    model.component('comp1').selection('box3').set('ymax', 'CellSize');
    model.component('comp1').selection('box3').set('zmin', 0);
    model.component('comp1').selection('box3').set('zmax', 'CellSize');
    model.component('comp1').selection('box3').set('condition', 'inside');
    
    model.component('comp1').selection.create('box4', 'Box');
    model.component('comp1').selection('box4').set('entitydim', 2);
    model.component('comp1').selection('box4').label('YN');
    model.component('comp1').selection('box4').set('xmin', 0);
    model.component('comp1').selection('box4').set('xmax', 'CellSize');
    model.component('comp1').selection('box4').set('ymin', 0);
    model.component('comp1').selection('box4').set('ymax', 0);
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
    model.component('comp1').selection('box5').set('zmin', 'CellSize');
    model.component('comp1').selection('box5').set('zmax', 'CellSize');
    model.component('comp1').selection('box5').set('condition', 'inside');
    
    model.component('comp1').selection.create('box6', 'Box');
    model.component('comp1').selection('box6').set('entitydim', 2);
    model.component('comp1').selection('box6').label('ZN');
    model.component('comp1').selection('box6').set('xmin', 0);
    model.component('comp1').selection('box6').set('xmax', 'CellSize');
    model.component('comp1').selection('box6').set('ymin', 0);
    model.component('comp1').selection('box6').set('ymax', 'CellSize');
    model.component('comp1').selection('box6').set('zmin', 0);
    model.component('comp1').selection('box6').set('zmax', 0);
    model.component('comp1').selection('box6').set('condition', 'inside');    
    
    model.component('comp1').physics('shell').prop('ShapeProperty').set('order_displacement', 1);
    
    model.component('comp1').physics('shell').feature('to1').set('d', [num2str(thickness),'[mm]']);

    model.component('comp1').physics('shell').prop('ShapeProperty').set('order_displacement', 1);

    model.component('comp1').physics('shell').create('disp1', 'Displacement2', 2);
    model.component('comp1').physics('shell').feature('disp1').selection.named('box1');
    model.component('comp1').physics('shell').feature('disp1').set('Direction', [1; 1; 1]);
    model.component('comp1').physics('shell').feature('disp1').set('U0', {'u_disp'; '0'; '0'});
    
    model.component('comp1').physics('shell').create('disp2', 'Displacement2', 2);
    model.component('comp1').physics('shell').feature('disp2').selection.named('box2');
    model.component('comp1').physics('shell').feature('disp2').set('Direction', [1; 1; 1]);
    model.component('comp1').physics('shell').feature('disp2').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('shell').create('disp3', 'Displacement2', 2);
    model.component('comp1').physics('shell').feature('disp3').selection.named('box3');
    model.component('comp1').physics('shell').feature('disp3').set('Direction', [0; 1; 1]);
    model.component('comp1').physics('shell').feature('disp3').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('shell').create('disp4', 'Displacement2', 2);
    model.component('comp1').physics('shell').feature('disp4').selection.named('box4');
    model.component('comp1').physics('shell').feature('disp4').set('Direction', [0; 1; 1]);
    model.component('comp1').physics('shell').feature('disp4').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('shell').create('disp5', 'Displacement2', 2);
    model.component('comp1').physics('shell').feature('disp5').selection.named('box5');
    model.component('comp1').physics('shell').feature('disp5').set('Direction', [0; 1; 1]);
    model.component('comp1').physics('shell').feature('disp5').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('shell').create('disp6', 'Displacement2', 2);
    model.component('comp1').physics('shell').feature('disp6').selection.named('box6');
    model.component('comp1').physics('shell').feature('disp6').set('Direction', [0; 1; 1]);
    model.component('comp1').physics('shell').feature('disp6').set('U0', {'0'; '0'; '0'});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    model.component('comp1').physics('shell').create('disp7', 'Displacement2', 2);
    model.component('comp1').physics('shell').feature('disp7').selection.named('box1');
    model.component('comp1').physics('shell').feature('disp7').set('Direction', [1; 0; 1]);
    model.component('comp1').physics('shell').feature('disp7').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('shell').create('disp8', 'Displacement2', 2);
    model.component('comp1').physics('shell').feature('disp8').selection.named('box2');
    model.component('comp1').physics('shell').feature('disp8').set('Direction', [1; 0; 1]);
    model.component('comp1').physics('shell').feature('disp8').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('shell').create('disp9', 'Displacement2', 2);
    model.component('comp1').physics('shell').feature('disp9').selection.named('box3');
    model.component('comp1').physics('shell').feature('disp9').set('Direction', [1; 1; 1]);
    model.component('comp1').physics('shell').feature('disp9').set('U0', {'0'; 'u_disp'; '0'});
    
    model.component('comp1').physics('shell').create('disp10', 'Displacement2', 2);
    model.component('comp1').physics('shell').feature('disp10').selection.named('box4');
    model.component('comp1').physics('shell').feature('disp10').set('Direction', [1; 1; 1]);
    model.component('comp1').physics('shell').feature('disp10').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('shell').create('disp11', 'Displacement2', 2);
    model.component('comp1').physics('shell').feature('disp11').selection.named('box5');
    model.component('comp1').physics('shell').feature('disp11').set('Direction', [1; 0; 1]);
    model.component('comp1').physics('shell').feature('disp11').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('shell').create('disp12', 'Displacement2', 2);
    model.component('comp1').physics('shell').feature('disp12').selection.named('box6');
    model.component('comp1').physics('shell').feature('disp12').set('Direction', [1; 0; 1]);
    model.component('comp1').physics('shell').feature('disp12').set('U0', {'0'; '0'; '0'});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    model.component('comp1').physics('shell').create('disp13', 'Displacement2', 2);
    model.component('comp1').physics('shell').feature('disp13').selection.named('box1');
    model.component('comp1').physics('shell').feature('disp13').set('Direction', [1; 1; 0]);
    model.component('comp1').physics('shell').feature('disp13').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('shell').create('disp14', 'Displacement2', 2);
    model.component('comp1').physics('shell').feature('disp14').selection.named('box2');
    model.component('comp1').physics('shell').feature('disp14').set('Direction', [1; 1; 0]);
    model.component('comp1').physics('shell').feature('disp14').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('shell').create('disp15', 'Displacement2', 2);
    model.component('comp1').physics('shell').feature('disp15').selection.named('box3');
    model.component('comp1').physics('shell').feature('disp15').set('Direction', [1; 1; 0]);
    model.component('comp1').physics('shell').feature('disp15').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('shell').create('disp16', 'Displacement2', 2);
    model.component('comp1').physics('shell').feature('disp16').selection.named('box4');
    model.component('comp1').physics('shell').feature('disp16').set('Direction', [1; 1; 0]);
    model.component('comp1').physics('shell').feature('disp16').set('U0', {'0'; '0'; '0'});
    
    model.component('comp1').physics('shell').create('disp17', 'Displacement2', 2);
    model.component('comp1').physics('shell').feature('disp17').selection.named('box5');
    model.component('comp1').physics('shell').feature('disp17').set('Direction', [1; 1; 1]);
    model.component('comp1').physics('shell').feature('disp17').set('U0', {'0'; '0'; 'u_disp'});
    
    model.component('comp1').physics('shell').create('disp18', 'Displacement2', 2);
    model.component('comp1').physics('shell').feature('disp18').selection.named('box6');
    model.component('comp1').physics('shell').feature('disp18').set('Direction', [1; 1; 1]);
    model.component('comp1').physics('shell').feature('disp18').set('U0', {'0'; '0'; '0'});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % model.save([Directory_files,name_of_model])
    
    % disp('Mesh')
    model.component('comp1').mesh.create('mesh1');
    model.component('comp1').mesh('mesh1').run;
    
    % model.save([Directory_files,name_of_model])

    % disp('Run')
    model.study.create('std1');
    model.study('std1').create('stat', 'Stationary');
    model.study('std1').feature('stat').set('useadvanceddisable', true);
    model.study('std1').feature('stat').set('disabledphysics', {'shell/disp7' 'shell/disp8' 'shell/disp9' 'shell/disp10' 'shell/disp11' 'shell/disp12' 'shell/disp13' 'shell/disp14' 'shell/disp15' 'shell/disp16'  ...
    'shell/disp17' 'shell/disp18'});
    
    model.study.create('std2');
    model.study('std2').create('stat', 'Stationary');
    model.study('std2').feature('stat').set('useadvanceddisable', true);
    model.study('std2').feature('stat').set('disabledphysics', {'shell/disp1' 'shell/disp2' 'shell/disp3' 'shell/disp4' 'shell/disp5' 'shell/disp6' 'shell/disp13' 'shell/disp14' 'shell/disp15' 'shell/disp16'  ...
    'shell/disp17' 'shell/disp18'});
    
    model.study.create('std3');
    model.study('std3').create('stat', 'Stationary');
    model.study('std3').feature('stat').set('useadvanceddisable', true);
    model.study('std3').feature('stat').set('disabledphysics', {'shell/disp1' 'shell/disp2' 'shell/disp3' 'shell/disp4' 'shell/disp5' 'shell/disp6' 'shell/disp7' 'shell/disp8' 'shell/disp9' 'shell/disp10'  ...
    'shell/disp11' 'shell/disp12'});
    
    model.save([Directory_files,name_of_model])
    
    model.study('std1').run;
    model.study('std2').run;
    model.study('std3').run;


    % model.save([Directory_files,name_of_model])

    UEnergy1=mphglobal(model,'shell.Ws_tot', 'dataset', 'dset1');
    UEnergy2=mphglobal(model,'shell.Ws_tot', 'dataset', 'dset2');
    UEnergy3=mphglobal(model,'shell.Ws_tot', 'dataset', 'dset3');
    
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