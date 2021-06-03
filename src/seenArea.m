function totalArea = seenArea(stlModel, seenMap)

seenFaces = stlModel.faces(logical(seenMap),:);

totalArea = 0;

    for ii = 1 : size(seenFaces)
        p1idx = seenFaces(ii, 1);
        p2idx = seenFaces(ii, 2);
        p3idx = seenFaces(ii, 3);
        vertices = [stlModel.vertices(p1idx, :); 
                    stlModel.vertices(p2idx, :);
                    stlModel.vertices(p3idx, :)];
        totalArea = totalArea + triangleArea(vertices);
    end
    
end