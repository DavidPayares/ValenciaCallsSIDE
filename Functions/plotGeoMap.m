function plotGeoMap(boundary, multipolygon)
    if (multipolygon == 1)
        x_min = findBoundingBoxValue(boundary ,1, 1) - 0.01;
        x_max = findBoundingBoxValue(boundary ,2, 0) + 0.01;
        y_min = findBoundingBoxValue(boundary ,3, 1) - 0.01;
        y_max = findBoundingBoxValue(boundary ,4, 0) + 0.01;
        ax = worldmap([y_min,y_max], [x_min,x_max]);
    else
        low = boundary.BoundingBox(1,:) - 0.01;
        high = boundary.BoundingBox(2,:) + 0.01;
        ax = worldmap([low(1,2),high(1,2)], [low(1,1),high(1,1)]);
    end
    setm(ax,'Grid','on','Frame','on','meridianlabel','on','parallellabel','on','FontSize',14)
    set(gcf,'color',[1,1,1])
    geoshow(ax, boundary, 'FaceColor', [0.93 0.87 0.51], 'FaceAlpha', 0.5) 
end