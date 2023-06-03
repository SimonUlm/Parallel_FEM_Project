clc
clear
close all

%% Define Constants
opt_         = struct('name',{},'val', {});
opt_(1).name = 'Node Number';
opt_(1).val  = 0;
opt_(2).name = 'Element Number';
opt_(2).val  = 0;
opt_(3).name = 'Edge Number';
opt_(3).val  = 0;
opt_(4).name = 'Affiliation';
opt_(4).val  = 0;
opt_(5).name = 'ComBorders';
opt_(5).val  = 0;
opt_(6).name = 'Boundary Elements';
opt_(6).val  = 0;
opt_(7).name = 'Fixed Notes';
opt_(7).val  = 0;
opt_(8).name = 'Print Elements';
opt_(8).val  = 0;

plot_only = 0;

opt = -1;
%% Debug-Mode
debugmode = 0; % 0 for user run --- 1 for debugging
opt_in_debugm = 1;
file_name = 'scattered_output_refinement_1.txt';


%% Select File for Mesh-Ploting
while(opt ~= 0)
    opt_file = -1;

    listing = dir('outputs');

    while(debugmode == 0 && opt_file ~= 0 && (all(opt_file ~= 3:length(listing))))
        clc
        for i = 1:length(listing)
            fprintf("["+i+"] "+listing(i).name+"\n")
        end
        opt_file = input(" ");

        if(isempty(opt_file))
            opt_file = -1;
        end
    end

    if opt_file > 0
        file_name = listing(opt_file).name;
    end

    if opt_file == 0
        break
    end

    %% read in .txt file
    file_text = readlines("outputs/" + file_name);
    file_text = strjoin(file_text);
    opt        = -1;   % User input for which thing to plot the numbers forz


    file_c_numMesh = length(strfind(file_text,'=========== Print Mesh Data ==========='));
    file_c_numSkel = length(strfind(file_text,'=========== Print Skeleton Data ==========='));
    
    file_text_mesh_array = extractBetween(file_text,...
    "=========== Print Mesh Data ===========",...
    "Memory");
    
    if(file_c_numSkel)
        file_text_Skel_array = extractBetween(file_text,...
        "=========== Print Skeleton Data ===========",...
        "Memory");
    end

    file_c_ComB = cell(file_c_numMesh,1); 
    file_c_Bele = cell(file_c_numMesh,1); 
    file_c_fixd = cell(file_c_numMesh,1); 
    M_coor      = cell(file_c_numMesh,1);
    M_elem      = cell(file_c_numMesh,1);
    M_Bord      = cell(file_c_numMesh,1);
    M_Bele      = cell(file_c_numMesh,1);
    M_BorN      = cell(file_c_numMesh,1);
    M_fixNod    = cell(file_c_numMesh,1);

    if(file_c_numSkel)
        for i = 1:file_c_numMesh
            [file_c_ComB{i}, file_c_Bele{i}, file_c_fixd{i}, M_coor{i},...
                M_elem{i}, M_Bord{i}, M_Bele{i}, M_BorN{i}, M_fixNod{i}] ...
                = text_extraction_2(file_text_mesh_array{i} + "Memory"...
                                  + file_text_Skel_array{i} + "Memory");
        end
    else
        for i = 1:file_c_numMesh
            [file_c_ComB{i}, file_c_Bele{i}, file_c_fixd{i}, M_coor{i},...
                M_elem{i}, M_Bord{i}, M_Bele{i}, M_BorN{i}, M_fixNod{i}] ...
                = text_extraction_2(file_text_mesh_array{i} + "Memory");
        end
    end

    %% Main Programm
    while (opt ~= 0)

        %% Reeding in user-input for which thing to plot
        % reset all variables
        opt = -1;
        vector_options = zeros(1, 9);
        vector_options(1:4) = (1:4);
        helper_for_opt = all(opt ~= vector_options);
        for i = 1:length(opt_)
            opt_(i).val = 0;
        end


        while (helper_for_opt && debugmode == 0)
            clc
            s_file_c = "";
            helper_for_opt = 0;


            if file_c_ComB{1}
                s_file_c = s_file_c + "\n[5] ComBorders";
                vector_options(5) = 5;
            end

            if file_c_Bele{1}
                s_file_c = s_file_c + "    [6] Boundary Elements";
                vector_options(6) = 6;
            end

            if file_c_fixd{1}
                s_file_c = s_file_c + "    [7] Fixed Notes ";
                vector_options(7) = 7;
            end

            if file_c_numMesh > 1
                s_file_c = s_file_c + "[8] which Element";
                vector_options(8) = 8;
            end
            file_c_numMesh
            disp(file_name)
            %%
            opt = input("[1] Node Number   [2] Element Number       " + ...
                "[3] Edge Number   [4] Affiliation"+s_file_c+"\n");

            if(isempty(opt))
                opt = -1;
            end

            for i = 1:length(opt)
                helper_for_opt = helper_for_opt + all(opt(i) ~= vector_options);
            end
        end

        if debugmode == 1
            opt = opt_in_debugm;
        end

        if (length(opt)== 1 && opt == 0)
            break
        end

        for i = 1:length(opt)
            opt_(opt(i)).val = 1;
        end
        
        if opt_(8).val
            plot_only = input("[0] all elements [n+1] n'th element\n");
        else
            if(~plot_only)
                for i = 1:file_c_numMesh
                    mesh_plot(file_c_ComB{i}, file_c_Bele{i}, file_c_fixd{i}, ...
                        M_coor{i}, M_elem{i}, M_Bord{i}, M_Bele{i},...
                        M_BorN{i}, M_fixNod{i}, opt_, true)
                end
            else
                for i = 1:file_c_numMesh
                    mesh_plot(file_c_ComB{i}, file_c_Bele{i}, file_c_fixd{i}, ...
                        M_coor{i}, M_elem{i}, M_Bord{i}, M_Bele{i},...
                        M_BorN{i}, M_fixNod{i}, opt_, any(i == plot_only))
                end
            end
            hold off
        end
        

        if(debugmode == 1)
            break;
        end

    end
    opt = -1;
end


