clc
clear
close all

opt = 10;
while (opt ~= 0)
    opt = 10;
    while (opt ~= 1 && opt ~= 2 && opt ~= 3 && opt ~= 4 && opt ~= 0)
        clc
        opt = input("[1] Node Number, [2] Element Number, [3] Edge Number, [4] Affiliation\n");
        
        if(isempty(opt))
            opt = 10;
        end
    end
    
    if (opt == 0)
        break
    end
    %% read in .txt file
    file_text = fileread('test_output.txt');

    Coordinates       = extractBetween(file_text,...
        "Coordinates (x,y):",...
        "Elements:");
    Elements          = extractBetween(file_text,...
        "Vertices (n1,n2,n3), Mid. Points (m1,m2,m3), Affiliation",...
        "Boundary Elements:");
    Boundary_Elements = extractBetween(file_text,...
        "Endpoints (n1, n2), Edge Number (ed1), Type",...
        "Fixed Nodes:");
    Fixed_Nodes       = extractBetween(file_text,"Fixed Nodes:", "Memory");
    %Fixed_Edges       = extractBetween(file_text,"Fixed Edges:", "Memory");

    Coordinates{1} = erase(Coordinates{1},["(", ",",")"]);
%     Fixed_Edges{1} = erase(Fixed_Edges{1},...
%         "Endpoints (n1, n2), Edge Number (ed1), Type");

    M_coor   = str2num(Coordinates{1});
    M_elem   = str2num(Elements{1});
    M_Bele   = str2num(Boundary_Elements{1});
    M_fixNod = str2num(Fixed_Nodes{1});
%     M_fixEdg = str2num(Fixed_Edges{1});


    %%
    plot(M_coor(:,1), M_coor(:,2), "o")
    hold on
    xlim([-0.5,1.5])
    ylim([-0.5,1.5])
    axis equal

    % print Node Number
    if(opt == 1)
        text(M_coor(:,1), M_coor(:,2), string(1:numel(M_coor(:,1))));
    end

    for j = 1:size(M_elem, 1)
        x_elem = zeros(1,3);
        y_elem = zeros(1,3);
        for i = 1:3
            x_elem(i) = M_coor(M_elem(j,i)+1,1);
            y_elem(i) = M_coor(M_elem(j,i)+1,2);
        end
        p = polyshape(x_elem, y_elem);
        plot(p);

        % print Element-Number
        if(opt == 2)
            text(mean(x_elem), mean(y_elem), string(j-1));
        end


        %print Edgenumber
        if(opt == 3)
            text(mean(x_elem(1:2)), mean(y_elem(1:2)), string(M_elem(j,4)));
            text(mean(x_elem(2:3)), mean(y_elem(2:3)), string(M_elem(j,5)));
            text(mean([x_elem(1),x_elem(3)]), mean([y_elem(1),y_elem(3)]), string(M_elem(j,6)));
        end

        %print Element-Affiliation
        if(opt == 4)
            text(mean(x_elem), mean(y_elem), string(M_elem(j,7)));
        end
        
        %print Boundary Elements
        

    end
    hold off
end


