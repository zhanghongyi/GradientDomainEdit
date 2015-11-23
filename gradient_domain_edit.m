function output = gradient_domain_edit(source, target, mask, loc)
    % process mask
    prop = regionprops(mask,'centroid', 'BoundingBox');
    % find center and rectangle box for mask
    center = cat(1, round(prop.Centroid));
    boundingbox = cat(1, round(prop.BoundingBox));
    % calculate mask coordinates
    y1 = center(1) - boundingbox(1) + 1;
    y2 = boundingbox(3) - y1 + 1;
    x1 = center(2) - boundingbox(2) + 1;
    x2 = boundingbox(4) - x1 + 1;
    % create source_mask and cropped masks
    cropped_mask = mask(center(2) - x1:center(2) + x2, center(1) - y1:center(1) + y2);
    cropped_target = target(center(2) - x1:center(2) + x2, center(1) - y1:center(1) + y2, :);

    s_size = size(cropped_mask);
    matrix_sz = s_size(1)*s_size(2);

    % construct sparse matrix
    % each pixel has 5 non-zero coeffiencts at most
    i=zeros(5*matrix_sz, 1);
    i(1:matrix_sz) = 1:matrix_sz;
    j=zeros(5*matrix_sz, 1);
    j(1:matrix_sz) = 1:matrix_sz;
    v=ones(5*matrix_sz, 1); %first matrix_sz number is reserved

    count = matrix_sz;
    for index=1:matrix_sz
        if(cropped_mask(index) == 1) %inside or boundary
            v(index) = 4;
            [x, y] = ind2sub([s_size(1), s_size(2)], index);  
            if(x > 1 && cropped_mask(x-1, y) == 1)
                count = count + 1;
                i(count) = index;
                j(count) = sub2ind([s_size(1), s_size(2)], x-1, y);
                v(count) = -1; 
            end
            if( y>1 && cropped_mask(x, y-1) == 1)
                count = count + 1;
                i(count) = index;
                j(count) = sub2ind([s_size(1), s_size(2)], x, y-1);
                v(count) = -1;    
            end
            if(x<size(cropped_mask, 1) && cropped_mask(x+1, y) == 1)
                count = count + 1;
                i(count) = index;
                j(count) = sub2ind([s_size(1), s_size(2)], x+1, y);
                v(count) = -1;
            end
            if(y < size(cropped_mask, 2) && cropped_mask(x, y+1) == 1)
                count = count + 1;
                i(count) = index;
                j(count) = sub2ind([s_size(1), s_size(2)], x, y+1);
                v(count) = -1;
            end
        end
    end
    i = i(1:count);
    j = j(1:count);
    v = v(1:count);
    A = sparse(i, j, v, matrix_sz, matrix_sz); % A is the same for every channel

    h = [ 0, -1, 0; -1, 4, -1; 0, -1, 0]; % laplacian matrix

    for n=1:3
        % compute for each channel
        b=ones(matrix_sz, 1);
        source_n = source(:, :, n);
        target_n = cropped_target(:, :, n);

        gradient_target_n = conv2(target_n, h, 'same');

        for index=1:matrix_sz
            if(cropped_mask(index) == 1) %inside or boundary
                % change value of right hand of the equation if the neighbor is
                % not in the mask
                b(index) = gradient_target_n(index);
                [x, y] = ind2sub([s_size(1), s_size(2)], index);  
                source_x = x - 1 + loc(1) - x1;
                source_y = y - 1 + loc(2) - y1;
                if(x == 1 || cropped_mask(x-1, y) == 0)
                    b(index) = b(index) + source_n(source_x-1, source_y);   
                end
                if( y==1 || cropped_mask(x, y-1) == 0)
                    b(index) = b(index) + source_n(source_x, source_y-1);     
                end
                if(x == size(target_n, 1) || cropped_mask(x+1, y) == 0)
                    b(index) = b(index) + source_n(source_x+1, source_y);  
                end
                if(y == size(target_n, 2) || cropped_mask(x, y+1) == 0)
                    b(index) = b(index) + source_n(source_x, source_y+1);   
                end
            end
        end
        
        % solve sparse matrix
        f = A\b;
        cropped_target(:, :, n) = reshape(f, size(target_n));
    end

    % full RGB
    output = source;
    output(loc(1) - x1:loc(1) + x2, loc(2) - y1:loc(2) + y2, :) = ...
        output(loc(1) - x1:loc(1) + x2, loc(2) - y1:loc(2) + y2, :) .* repmat(1-cropped_mask, 1, 1, 3) + ...
        cropped_target .* repmat(cropped_mask, 1, 1, 3);
end