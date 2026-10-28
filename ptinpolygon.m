function inside_points = ptinpolygon(x,y,pointx,pointy)
polygons = [x, y];
is_inside = false(length(pointx(:)),1);
for k = 1:length(pointx)
    xq = pointx(k);
    yq = pointy(k);

    [in, ~] = inpolygon(xq, yq, polygons(:,1), polygons(:,2));
    if in
        is_inside(k,1) = true;
    end
end
inside_points = [pointx(is_inside), pointy(is_inside)];