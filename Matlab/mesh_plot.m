function mesh_plot(file_c_ComB, file_c_Bele, file_c_fixd, M_coor, M_elem, M_Bord, M_Bele, M_BorN, M_fixNod, opt_, just_plot)
color      = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'w'];
color4     = ['r', 'b', 'k', 'm'];
l_style    = {'--', ':'};



%% Plotting All Nodes
plot(M_coor(:,1), M_coor(:,2), ".")
hold on
xlim([-0.5,1.5])
ylim([-0.5,1.5])
axis equal

%% print Node Numbers
if(opt_(1).val == 1 && just_plot)
    text(M_coor(:,1), M_coor(:,2), string(0:numel(M_coor(:,1))-1));
end

%% Plotting All Colored ComBorders

if((opt_(5).val == 1) && file_c_ComB)
    for i = 1:size(M_Bord,1)
        my_color = color4(M_Bord(i,6)+1);

        plot(M_coor(M_Bord(i,2)+1,1), M_coor(M_Bord(i,2)+1,2),...
            'd', 'MarkerSize', 10, 'MarkerEdgeColor', my_color,...
            'MarkerFaceColor', my_color)
        plot(M_coor(M_Bord(i,3)+1,1), M_coor(M_Bord(i,3)+1,2),...
            'd', 'MarkerSize', 10, 'MarkerEdgeColor', my_color,...
            'MarkerFaceColor', my_color)

        x = M_coor(M_BorN(M_Bord(i,1)+1,:)'+1, 1);
        y = M_coor(M_BorN(M_Bord(i,1)+1,:)'+1, 2);

        plot(x, y,...
            's', 'MarkerSize', 10, 'MarkerEdgeColor', my_color,...
            'MarkerFaceColor', my_color)
    end
end


%% Plotting regarding Elements
for j = 1:size(M_elem, 1)
    %Plot Elements
    x_elem = zeros(1,3);
    y_elem = zeros(1,3);
    for i = 1:3
        x_elem(i) = M_coor(M_elem(j,i)+1,1);
        y_elem(i) = M_coor(M_elem(j,i)+1,2);
    end
    p = polyshape(x_elem, y_elem);
    plot(p,'FaceColor', color(mod(M_elem(j,7), 8)+1));

    % print Element-Number
    if(opt_(2).val == 1)
        text(mean(x_elem), mean(y_elem), string(j-1));
    end


    % print Edgenumber
    if(opt_(3).val == 1 && just_plot)
        text(mean(x_elem(1:2)), mean(y_elem(1:2)), string(M_elem(j,4)));
        text(mean(x_elem(2:3)), mean(y_elem(2:3)), string(M_elem(j,5)));
        text(mean([x_elem(1),x_elem(3)]), mean([y_elem(1),y_elem(3)]), string(M_elem(j,6)));
    end

    %print Element-Affiliation
    if(opt_(4).val == 1)
        text(mean(x_elem), mean(y_elem), string(M_elem(j,7)));
    end
end

%print Boundary Elements
if((opt_(6).val == 1) && (file_c_Bele))
    for i = 1:size(M_Bele,1)
        x = [M_coor(M_Bele(i,1)+1,1), M_coor(M_Bele(i,2)+1,1)];
        y = [M_coor(M_Bele(i,1)+1,2), M_coor(M_Bele(i,2)+1,2)];
        plot(x, y, l_style{M_Bele(i,4)+1}, 'LineWidth',3, ...
            'Color', 	"#7E2F8E")
        text(mean(x), ...
            mean(y), ...
            string(M_Bele(i,3)));
    end
end



%print fixed Notes
if((opt_(7).val == 1) && file_c_fixd)
    for i = 1:size(M_fixNod,1)
        plot(M_coor(M_fixNod(i,1)+1,1), M_coor(M_fixNod(i,1)+1,2),...
            'o', 'MarkerSize', 10, 'MarkerEdgeColor', "#7E2F8E",...
            'MarkerFaceColor', 	"#7E2F8E")
        text(M_coor(M_fixNod(i,1)+1,1), M_coor(M_fixNod(i,1)+1,2), string(M_fixNod(i,1)));
    end
end
end