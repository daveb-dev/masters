% for i = 2:40
%     for j = 2:60
%         for k = 2:15
%             if skull(i,j,k) == 1
%                 count = 0;
%                 for a = i-1:i+1
%                     for b = j-1:j+1
%                         for c = k-1:k+1
%                             if skull(a,b,c) == 1
%                                 count = count + 1;
%                             end
%                         end
%                     end
%                 end
%                 if count <= 2
%                     fprintf('%d,%d,%d',i,j,k);
%                     break
%                 else
%                     fprintf('%d,%d,%d,%i\n',i,j,k,count);
%                 end
%             end
%         end
%     end
% end
close

x = []; y = []; z = 7;
for i = 1:41
    for j = 1:61
        if skull(i,j,z) > 0
            x = [x; i]; y = [y; j];
        end
    end
end
scatter(x,y); hold on

x = []; y = []; c = []; t = 3;
for i = 1:41
    for j = 1:61
        if cells(i,j,z,t) > 0
            x = [x; i]; y = [y; j]; c = [c; cells(i,j,z,t)];
        end
    end
end
scatter(x,y,[],c)