function val = findBoundingBoxValue(shapefile ,corner, minimum)

vals = [];
for i = 1:length(shapefile)
    num = shapefile(i).BoundingBox(corner);
    vals = [vals, num];
end

if minimum == 1
    val = min(vals);
else
    val = max(vals);
end

end