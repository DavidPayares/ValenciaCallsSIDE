function vp = plotCovariateMap(boundary, x_field, y_field, covariate_name, mask, S1, S2 , S1_lat, S2_lon)
    vp = griddata(x_field, y_field, covariate_name, S1,S2);
    figure()
    plotGeoMap(boundary,0)
    vp = vp.*mask;
    surfm(S1_lat,S2_lon,vp)
    shading interp
    c = flipud(hot);
    colormap(c)
    colorbar
    caxis([0, max(max(vp))]);
end