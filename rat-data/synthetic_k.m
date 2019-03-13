clear

cd(strcat('./rat05'));    % Load skull boundary and initial tumor
load('skull_out.mat');
load('tumor_out.mat');
load('tumor_t0.mat');
cd ..

[m,n] = size(tumor); % Known size of domain
k = zeros(m,n);
tumor_out(:,2) = 41-tumor_out(:,2); % This was flipped in plotdata
skull_out(:,2) = 41-skull_out(:,2);

tumor_min_x = min(tumor_out(:,2)); % Tumor bounds
tumor_max_x = max(tumor_out(:,2));
tumor_min_y = min(tumor_out(:,1));
tumor_max_y = max(tumor_out(:,1));

max_distance = 6;

min_rand_x = max(1,tumor_min_x - max_distance);
max_rand_x = min(m,tumor_max_x + max_distance);
min_rand_y = max(1,tumor_min_y - max_distance);
max_rand_y = min(n,tumor_max_y + max_distance);

% Select coordinates for random source location
x = min_rand_x+round((max_rand_x-min_rand_x)*rand());     
y = min_rand_y+round((max_rand_y-min_rand_y)*rand());   

max_radius = 9;

for i = max_radius:-1:1
    k_value = 2*(max_radius-i)/max_radius;
    x_lower_bound = max(1,x-i);
    x_upper_bound = min(m,x+i);
    y_lower_bound = max(1,y-i);
    y_upper_bound = min(n,y+i);
    for j = x_lower_bound:x_upper_bound
        for l = y_lower_bound:y_upper_bound
            a = abs(x-j);
            b = abs(y-l);
            c = sqrt(a^2+b^2);
            if c < i
                k(j,l) = k_value;
            end
        end
    end
end
k(x,y) = 2;

I = find(k>0); % Uses linear indexing to find coordinates
xq = mod(I,41);
yq = ceil(I/41);
[in] = inpolygon(xq, yq, skull_out(:,2), skull_out(:,1));

xq = xq(~in);
yq = yq(~in);
for i = 1:length(xq)
    k(xq(i),yq(i)) = 0;
end

imshow(k)
save("./k_synth.mat",'k');