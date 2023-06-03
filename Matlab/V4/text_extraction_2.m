function [file_c_ComB, file_c_Bele, file_c_fixd, M_coor, M_elem, M_Bord, M_Bele, M_BorN, M_fixNod] = text_extraction_2(file_text)

%% Check file for Keywords
file_c_ComB  = contains(file_text, 'ComBorders:');  % File contains Color opt
file_c_Bele  = contains(file_text, 'Boundary Elements:');
file_c_fixd  = contains(file_text, 'Fixed Nodes:');

%% Determin the end of the specific sections
if(file_c_fixd)
    End_boundary_Elements = "Fixed Nodes:";
else
    End_boundary_Elements =  "Memory";
end

%% Extract data from .txt file

Coordinates       = extractBetween(file_text,...
    "Coordinates (x,y):",...
    "Elements:");
Coordinates{1} = erase(Coordinates{1},["(", ",",")"]);

%"Elements:" + newline + ...
Elements          = extractBetween(file_text,...
    "Vertices (n1,n2,n3), Mid. Points (m1,m2,m3), Affiliation",...
    "Boundary Elements:");
%"Boundary Elements:" + newline + ...
Boundary_Elements = extractBetween(file_text,...
    "Endpoints (n1, n2), Edge Number (ed1), Type",...
    End_boundary_Elements);

if( file_c_fixd )
    Fixed_Nodes       = extractBetween(file_text,"Fixed Nodes:", "Memory");
end

if( file_c_ComB == 1)
        %"ComBorders:" + newline + ...
    Borders          = extractBetween(file_text,...
        "Index, Start_Node, End_Node, L-Proc, R-Proc, Color",...
        "ComBorderNodes:");
    
        %"ComBorderNodes:" + newline + ...
    Border_Nodes     = extractBetween(file_text,...
        "Connecting Nodes",...
        "Memory");
end




M_coor   = str2num(Coordinates{1});
M_coor   = reshape(M_coor, 2,[])';
M_elem   = str2num(Elements{1});
M_elem   = reshape(M_elem, 7,[])';

if(file_c_ComB)
    M_Bord   = str2num(Borders{1});
    M_Bord   = reshape(M_Bord, 6,[])';

    M_BorN   = str2num(Border_Nodes{1});
    M_BorN   = reshape(M_BorN, [], size(M_Bord, 1))';
else
    M_Bord   = 'not in use';
    M_BorN   = 'not in use';
end

if(file_c_Bele)
    M_Bele   = str2num(Boundary_Elements{1});
    M_Bele   = reshape(M_Bele, 4,[])';
else
    M_Bele   = 'not in use';
end

if(file_c_fixd)
    M_fixNod = str2num(Fixed_Nodes{1});
    M_fixNod   = reshape(M_fixNod, 1,[])';
else
    M_fixNod = 'not in use';
end