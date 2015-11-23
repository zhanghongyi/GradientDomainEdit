function output = copy_paste_edit(source, target, mask, loc)
    % directly copy paste
    prop = regionprops(mask,'centroid', 'BoundingBox');
    center = cat(1, round(prop.Centroid));
    boundingbox = cat(1, round(prop.BoundingBox));
    y1 = center(1) - boundingbox(1) + 1;
    y2 = boundingbox(3) - y1 + 1;
    x1 = center(2) - boundingbox(2) + 1;
    x2 = boundingbox(4) - x1 + 1;
    source_mask = ones(size(source, 1), size(source, 2));
    cropped_mask = mask(center(2) - x1:center(2) + x2, center(1) - y1:center(1) + y2);
    cropped_target = target(center(2) - x1:center(2) + x2, center(1) - y1:center(1) + y2, :);

    source_mask(loc(1) - x1:loc(1) + x2, loc(2) - y1:loc(2) + y2) = 1 - cropped_mask;
    output = source .* repmat(source_mask, 1,1,3);
    output(loc(1) - x1:loc(1) + x2, loc(2) - y1:loc(2) + y2, :) = ...
        output(loc(1) - x1:loc(1) + x2, loc(2) - y1:loc(2) + y2, :) + ...
        cropped_target .* repmat(cropped_mask, 1, 1, 3);
end