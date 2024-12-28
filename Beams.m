function Beams(rownumber, columnnumber, cube_size, u_disp, rho, E, nu)

    % Directory_files='/home/rkm41/GraphDataZ/';
    Directory_files='C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\';
    NameData = ['New_outer_Beam',num2str(rownumber)];
    import com.comsol.model.*
    import com.comsol.model.util.*

    StartNodes_all = [];
    EndNodes_all = [];
    Number_of_Nodes=0;

    name_of_model=[NameData,'.mph'];
    for i=1:rownumber
        for j=1:columnnumber
            NameData = ['Rod',num2str((i-1)*rownumber+j)];
            nodes=load([Directory_files,NameData,'.mat']);
            Number_of_Node=length(nodes.edges);
            StartNodes=zeros(Number_of_Node,3);
            StartNodes(:,:)=nodes.edges(:,1,:);
            EndNodes=zeros(Number_of_Node,3);
            EndNodes(:,:)=nodes.edges(:,2,:);

            Number_of_Nodes=Number_of_Nodes+Number_of_Node;

            StartNodes(:,1) = StartNodes(:,1) + cube_size*(j-1);
            StartNodes(:,2) = StartNodes(:,2) + cube_size*(i-1);
            EndNodes(:,1) = EndNodes(:,1) + cube_size*(j-1);
            EndNodes(:,2) = EndNodes(:,2) + cube_size*(i-1);

            StartNodes_all = [StartNodes_all; StartNodes(:,:)];
            EndNodes_all = [EndNodes_all; EndNodes(:,:)];
        end
    end
    
    [~,I,~]=unique(StartNodes_all,'rows');
    Index_Identique_Start=setdiff(1:length(StartNodes_all),I);
    
    [~,J,~]=unique(EndNodes_all,'rows');
    Index_Identique_End=setdiff(1:length(EndNodes_all),J);
    
    r_cylinder=0.01; %mm
    
    axis_cylinder = EndNodes_all - StartNodes_all;
    l_cylinder = vecnorm(axis_cylinder, 2, 2);
    
    model = ModelUtil.create('Model');
    
    model.param.set('CellSize', [num2str(cube_size), ' [mm]']);
    model.param.set('u_disp', [num2str(u_disp), ' [mm]']);
    model.param.set('r_cylinder', [num2str(r_cylinder), ' [mm]']);
    model.param.set('rownumber', num2str(rownumber));
    model.param.set('columnnumber', num2str(columnnumber));
    
    model.component.create('comp1', true);
    
    % disp('Geometry')
    
    model.component('comp1').geom.create('geom1', 3);
    
    model.component('comp1').geom('geom1').lengthUnit('mm');
    % model.component('comp1').geom('geom1').geomRep('comsol');
    model.component('comp1').geom('geom1').geomRep('cadps');


    for i=1:Number_of_Nodes
        % disp(i)
        model.component('comp1').geom('geom1').create(['ls',num2str(i)], 'LineSegment');
        model.component('comp1').geom('geom1').feature(['ls',num2str(i)]).set('specify1', 'coord');
        model.component('comp1').geom('geom1').feature(['ls',num2str(i)]).set('coord1', [StartNodes_all(i,1) StartNodes_all(i,2) StartNodes_all(i,3)]);
        model.component('comp1').geom('geom1').feature(['ls',num2str(i)]).set('specify2', 'coord');
        model.component('comp1').geom('geom1').feature(['ls',num2str(i)]).set('coord2', [EndNodes_all(i,1) EndNodes_all(i,2) EndNodes_all(i,3)]);
    end

    model.component('comp1').geom('geom1').run;
    
    % model.save([Directory_files,name_of_model])
    
    model.component('comp1').massProp.create('mass1', 'MassProperties');
    model.component('comp1').massProp('mass1').set('densitySource', 'fromSpecifiedPhysics');
    % 
    model.component('comp1').selection.create('box1', 'Box');
    model.component('comp1').selection('box1').set('entitydim', 0);
    model.component('comp1').selection('box1').label('XP');
    model.component('comp1').selection('box1').set('xmin', 'CellSize*rownumber');
    model.component('comp1').selection('box1').set('xmax', 'CellSize*rownumber');
    model.component('comp1').selection('box1').set('ymin', 0);
    model.component('comp1').selection('box1').set('ymax', 'CellSize*rownumber');
    model.component('comp1').selection('box1').set('zmin', 0);
    model.component('comp1').selection('box1').set('zmax', 'CellSize');
    model.component('comp1').selection('box1').set('condition', 'inside');
    
    model.component('comp1').selection.create('box2', 'Box');
    model.component('comp1').selection('box2').set('entitydim', 0);
    model.component('comp1').selection('box2').label('XN');
    model.component('comp1').selection('box2').set('xmin', 0);
    model.component('comp1').selection('box2').set('xmax', 0);
    model.component('comp1').selection('box2').set('ymin', 0);
    model.component('comp1').selection('box2').set('ymax', 'CellSize*rownumber');
    model.component('comp1').selection('box2').set('zmin', 0);
    model.component('comp1').selection('box2').set('zmax', 'CellSize');
    model.component('comp1').selection('box2').set('condition', 'inside');
    
    model.component('comp1').selection.create('box3', 'Box');
    model.component('comp1').selection('box3').set('entitydim', 0);
    model.component('comp1').selection('box3').label('YP');
    model.component('comp1').selection('box3').set('xmin', 0);
    model.component('comp1').selection('box3').set('xmax', 'CellSize*rownumber');
    model.component('comp1').selection('box3').set('ymin', 'CellSize*rownumber');
    model.component('comp1').selection('box3').set('ymax', 'CellSize*rownumber');
    model.component('comp1').selection('box3').set('zmin', 0);
    model.component('comp1').selection('box3').set('zmax', 'CellSize');
    model.component('comp1').selection('box3').set('condition', 'inside');
    
    model.component('comp1').selection.create('box4', 'Box');
    model.component('comp1').selection('box4').set('entitydim', 0);
    model.component('comp1').selection('box4').label('YN');
    model.component('comp1').selection('box4').set('xmin', 0);
    model.component('comp1').selection('box4').set('xmax', 'CellSize*rownumber');
    model.component('comp1').selection('box4').set('ymin', 0);
    model.component('comp1').selection('box4').set('ymax', 0);
    model.component('comp1').selection('box4').set('zmin', 0);
    model.component('comp1').selection('box4').set('zmax', 'CellSize');
    model.component('comp1').selection('box4').set('condition', 'inside');
    
    model.component('comp1').selection.create('box5', 'Box');
    model.component('comp1').selection('box5').set('entitydim', 0);
    model.component('comp1').selection('box5').label('ZP');
    model.component('comp1').selection('box5').set('xmin', 0);
    model.component('comp1').selection('box5').set('xmax', 'CellSize*rownumber');
    model.component('comp1').selection('box5').set('ymin', 0);
    model.component('comp1').selection('box5').set('ymax', 'CellSize*rownumber');
    model.component('comp1').selection('box5').set('zmin', 'CellSize');
    model.component('comp1').selection('box5').set('zmax', 'CellSize');
    model.component('comp1').selection('box5').set('condition', 'inside');
    
    model.component('comp1').selection.create('box6', 'Box');
    model.component('comp1').selection('box6').set('entitydim', 0);
    model.component('comp1').selection('box6').label('ZN');
    model.component('comp1').selection('box6').set('xmin', 0);
    model.component('comp1').selection('box6').set('xmax', 'CellSize*rownumber');
    model.component('comp1').selection('box6').set('ymin', 0);
    model.component('comp1').selection('box6').set('ymax', 'CellSize*rownumber');
    model.component('comp1').selection('box6').set('zmin', 0);
    model.component('comp1').selection('box6').set('zmax', 0);
    model.component('comp1').selection('box6').set('condition', 'inside');
    
    model.component('comp1').material.create('mat1', 'Common');
    model.component('comp1').material('mat1').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');

    model.component('comp1').physics.create('beam', 'HermitianBeam', 'geom1');
    model.component('comp1').physics('beam').create('pdr1', 'DispRot0', 0);
    model.component('comp1').physics('beam').feature('pdr1').selection.named('box1');
    model.component('comp1').physics('beam').feature('pdr1').setIndex('Direction', 'prescribed', 0);
    model.component('comp1').physics('beam').feature('pdr1').setIndex('U0', '0.1 [mm]', 0);
    model.component('comp1').physics('beam').feature('pdr1').setIndex('Direction', 'prescribed', 1);
    model.component('comp1').physics('beam').feature('pdr1').setIndex('Direction', 'prescribed', 2);
    model.component('comp1').physics('beam').feature.duplicate('pdr2', 'pdr1');
    model.component('comp1').physics('beam').feature('pdr2').setIndex('U0', 0, 0);
    model.component('comp1').physics('beam').feature('pdr2').selection.named('box2');

    model.component('comp1').material('mat1').label('PR48 Acrylic photoresist');
    model.component('comp1').material('mat1').propertyGroup('def').set('density', num2str(rho));
    model.component('comp1').material('mat1').propertyGroup('Enu').set('E', [num2str(E), ' [GPa]']);
    model.component('comp1').material('mat1').propertyGroup('Enu').set('nu', num2str(nu));
    model.save([Directory_files,name_of_model])

    model.component('comp1').physics('beam').feature('csd1').set('SectionType', 'CircularSection');
    model.component('comp1').physics('beam').feature('csd1').set('do_circ', 'r_cylinder');
    model.component('comp1').physics('beam').feature('csd1').feature('so1').set('point_beam', [1000; 5000; 7000]);


    model.component('comp1').massProp('mass1').set('densitySource', 'fromSpecifiedPhysics');
    
    % model.component('comp1').physics('beam').prop('AdvancedSettings').set('GroupPhysOdesRd', false);
    % model.component('comp1').physics('beam').feature('csd1').set('CrossSectionDefinition', 'CommonSections');
    % model.component('comp1').physics('beam').feature('csd1').set('SectionType', 'CircularSection');
    % model.component('comp1').physics('beam').feature('csd1').set('do_circ', 'radius');
    % model.component('comp1').physics('beam').feature('csd1').feature('so1').set('point_beam', [100; 500; 700]);
    % model.component('comp1').physics('beam').feature('pdr1').set('Direction', {'prescribed'; 'prescribed'; 'prescribed'});
    % model.component('comp1').physics('beam').feature('pdr2').set('Direction', {'prescribed'; 'prescribed'; 'prescribed'});
    % model.component('comp1').physics('beam').feature('pdr2').set('U0', [1; 0; 0]);

    % 
    % model.component('comp1').physics('solid').create('disp1', 'Displacement2', 2);
    % model.component('comp1').physics('solid').feature('disp1').selection.named('box1');
    % model.component('comp1').physics('solid').feature('disp1').set('Direction', [1; 1; 1]);
    % model.component('comp1').physics('solid').feature('disp1').set('U0', {'u_disp'; '0'; '0'});
    % 
    % model.component('comp1').physics('solid').create('disp2', 'Displacement2', 2);
    % model.component('comp1').physics('solid').feature('disp2').selection.named('box2');
    % model.component('comp1').physics('solid').feature('disp2').set('Direction', [1; 1; 1]);
    % model.component('comp1').physics('solid').feature('disp2').set('U0', {'0'; '0'; '0'});
    % 
    % model.component('comp1').physics('solid').create('disp3', 'Displacement2', 2);
    % model.component('comp1').physics('solid').feature('disp3').selection.named('box3');
    % model.component('comp1').physics('solid').feature('disp3').set('Direction', [0; 1; 1]);
    % model.component('comp1').physics('solid').feature('disp3').set('U0', {'0'; '0'; '0'});
    % 
    % model.component('comp1').physics('solid').create('disp4', 'Displacement2', 2);
    % model.component('comp1').physics('solid').feature('disp4').selection.named('box4');
    % model.component('comp1').physics('solid').feature('disp4').set('Direction', [0; 1; 1]);
    % model.component('comp1').physics('solid').feature('disp4').set('U0', {'0'; '0'; '0'});
    % 
    % model.component('comp1').physics('solid').create('disp5', 'Displacement2', 2);
    % model.component('comp1').physics('solid').feature('disp5').selection.named('box5');
    % model.component('comp1').physics('solid').feature('disp5').set('Direction', [0; 1; 1]);
    % model.component('comp1').physics('solid').feature('disp5').set('U0', {'0'; '0'; '0'});
    % 
    % model.component('comp1').physics('solid').create('disp6', 'Displacement2', 2);
    % model.component('comp1').physics('solid').feature('disp6').selection.named('box6');
    % model.component('comp1').physics('solid').feature('disp6').set('Direction', [0; 1; 1]);
    % model.component('comp1').physics('solid').feature('disp6').set('U0', {'0'; '0'; '0'});
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % model.component('comp1').physics('solid').create('disp7', 'Displacement2', 2);
    % model.component('comp1').physics('solid').feature('disp7').selection.named('box1');
    % model.component('comp1').physics('solid').feature('disp7').set('Direction', [1; 0; 1]);
    % model.component('comp1').physics('solid').feature('disp7').set('U0', {'0'; '0'; '0'});
    % 
    % model.component('comp1').physics('solid').create('disp8', 'Displacement2', 2);
    % model.component('comp1').physics('solid').feature('disp8').selection.named('box2');
    % model.component('comp1').physics('solid').feature('disp8').set('Direction', [1; 0; 1]);
    % model.component('comp1').physics('solid').feature('disp8').set('U0', {'0'; '0'; '0'});
    % 
    % model.component('comp1').physics('solid').create('disp9', 'Displacement2', 2);
    % model.component('comp1').physics('solid').feature('disp9').selection.named('box3');
    % model.component('comp1').physics('solid').feature('disp9').set('Direction', [1; 1; 1]);
    % model.component('comp1').physics('solid').feature('disp9').set('U0', {'0'; 'u_disp'; '0'});
    % 
    % model.component('comp1').physics('solid').create('disp10', 'Displacement2', 2);
    % model.component('comp1').physics('solid').feature('disp10').selection.named('box4');
    % model.component('comp1').physics('solid').feature('disp10').set('Direction', [1; 1; 1]);
    % model.component('comp1').physics('solid').feature('disp10').set('U0', {'0'; '0'; '0'});
    % 
    % model.component('comp1').physics('solid').create('disp11', 'Displacement2', 2);
    % model.component('comp1').physics('solid').feature('disp11').selection.named('box5');
    % model.component('comp1').physics('solid').feature('disp11').set('Direction', [1; 0; 1]);
    % model.component('comp1').physics('solid').feature('disp11').set('U0', {'0'; '0'; '0'});
    % 
    % model.component('comp1').physics('solid').create('disp12', 'Displacement2', 2);
    % model.component('comp1').physics('solid').feature('disp12').selection.named('box6');
    % model.component('comp1').physics('solid').feature('disp12').set('Direction', [1; 0; 1]);
    % model.component('comp1').physics('solid').feature('disp12').set('U0', {'0'; '0'; '0'});
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % model.component('comp1').physics('solid').create('disp13', 'Displacement2', 2);
    % model.component('comp1').physics('solid').feature('disp13').selection.named('box1');
    % model.component('comp1').physics('solid').feature('disp13').set('Direction', [1; 1; 0]);
    % model.component('comp1').physics('solid').feature('disp13').set('U0', {'0'; '0'; '0'});
    % 
    % model.component('comp1').physics('solid').create('disp14', 'Displacement2', 2);
    % model.component('comp1').physics('solid').feature('disp14').selection.named('box2');
    % model.component('comp1').physics('solid').feature('disp14').set('Direction', [1; 1; 0]);
    % model.component('comp1').physics('solid').feature('disp14').set('U0', {'0'; '0'; '0'});
    % 
    % model.component('comp1').physics('solid').create('disp15', 'Displacement2', 2);
    % model.component('comp1').physics('solid').feature('disp15').selection.named('box3');
    % model.component('comp1').physics('solid').feature('disp15').set('Direction', [1; 1; 0]);
    % model.component('comp1').physics('solid').feature('disp15').set('U0', {'0'; '0'; '0'});
    % 
    % model.component('comp1').physics('solid').create('disp16', 'Displacement2', 2);
    % model.component('comp1').physics('solid').feature('disp16').selection.named('box4');
    % model.component('comp1').physics('solid').feature('disp16').set('Direction', [1; 1; 0]);
    % model.component('comp1').physics('solid').feature('disp16').set('U0', {'0'; '0'; '0'});
    % 
    % model.component('comp1').physics('solid').create('disp17', 'Displacement2', 2);
    % model.component('comp1').physics('solid').feature('disp17').selection.named('box5');
    % model.component('comp1').physics('solid').feature('disp17').set('Direction', [1; 1; 1]);
    % model.component('comp1').physics('solid').feature('disp17').set('U0', {'0'; '0'; 'u_disp'});
    % 
    % model.component('comp1').physics('solid').create('disp18', 'Displacement2', 2);
    % model.component('comp1').physics('solid').feature('disp18').selection.named('box6');
    % model.component('comp1').physics('solid').feature('disp18').set('Direction', [1; 1; 1]);
    % model.component('comp1').physics('solid').feature('disp18').set('U0', {'0'; '0'; '0'});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % model.save([Directory_files,name_of_model])
    
    % disp('Mesh')
    model.component('comp1').mesh.create('mesh1');
    model.component('comp1').mesh('mesh1').run;
    
    model.save([Directory_files,name_of_model])
    % 
    % % disp('Run')
    % model.study.create('std1');
    % model.study('std1').create('stat', 'Stationary');
    % model.study('std1').feature('stat').set('useadvanceddisable', true);
    % model.study('std1').feature('stat').set('disabledphysics', {'solid/disp7' 'solid/disp8' 'solid/disp9' 'solid/disp10' 'solid/disp11' 'solid/disp12' 'solid/disp13' 'solid/disp14' 'solid/disp15' 'solid/disp16'  ...
    % 'solid/disp17' 'solid/disp18'});
    % 
    % model.study.create('std2');
    % model.study('std2').create('stat', 'Stationary');
    % model.study('std2').feature('stat').set('useadvanceddisable', true);
    % model.study('std2').feature('stat').set('disabledphysics', {'solid/disp1' 'solid/disp2' 'solid/disp3' 'solid/disp4' 'solid/disp5' 'solid/disp6' 'solid/disp13' 'solid/disp14' 'solid/disp15' 'solid/disp16'  ...
    % 'solid/disp17' 'solid/disp18'});
    % 
    % model.study.create('std3');
    % model.study('std3').create('stat', 'Stationary');
    % model.study('std3').feature('stat').set('useadvanceddisable', true);
    % model.study('std3').feature('stat').set('disabledphysics', {'solid/disp1' 'solid/disp2' 'solid/disp3' 'solid/disp4' 'solid/disp5' 'solid/disp6' 'solid/disp7' 'solid/disp8' 'solid/disp9' 'solid/disp10'  ...
    % 'solid/disp11' 'solid/disp12'});
    % 
    % model.save([Directory_files,name_of_model])
    % 
    % model.study('std1').run;
    % model.study('std2').run;
    % model.study('std3').run;
    % 
    % % model.save([Directory_files,name_of_model])
    % 
    % UEnergy1=mphglobal(model,'solid.Ws_tot', 'dataset', 'dset1');
    % UEnergy2=mphglobal(model,'solid.Ws_tot', 'dataset', 'dset2');
    % UEnergy3=mphglobal(model,'solid.Ws_tot', 'dataset', 'dset3');
    % 
    % % U_data2=mpheval(model,'solid.Ws_tot', 'dataset', 'dset2');
    % % UEnergy2=U_data2.d1(1);
    % 
    % C11=2*UEnergy1/(cube_size/1000)^3/(u_disp/cube_size)^2;
    % C22=2*UEnergy2/(cube_size/1000)^3/(u_disp/cube_size)^2;
    % C33=2*UEnergy3/(cube_size/1000)^3/(u_disp/cube_size)^2;
    % 
    % MassValue=mphglobal(model,'mass1.mass');
    % DensityStructure=MassValue/(cube_size/1000)^3;
    % 
    % Impedance1 = sqrt(C11*DensityStructure);
    % Impedance2 = sqrt(C22*DensityStructure);
    % Impedance3 = sqrt(C33*DensityStructure);
    % 
    % VolumeFraction = MassValue/(rho*(cube_size/1000)^3);
    % 
    % Stiffness = [C11, C22, C33];
    % Impedance = [Impedance1, Impedance2, Impedance3];
    % 
    % output = struct();
    % output.Impedance = Impedance;
    % output.Stiffness = Stiffness;
    % output.VolumeFraction = VolumeFraction;
    % output.DensityStructure = DensityStructure;
    % 
    % % Save the new parameters to the existing .mat file, appending the data
    % save([Directory_files, NameData, '.mat'], '-struct', 'output', '-append');
    % % disp([Directory_files, NameData, '.mat'])

end
