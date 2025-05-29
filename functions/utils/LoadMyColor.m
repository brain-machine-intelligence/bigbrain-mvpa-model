function color = LoadMyColor(N)

if nargin == 0
    N = 10;
end

% elseif nargin > 1
%         error('Check your input');
% else
%     if N > 10 || N < 1
%         N = 10;
%     end
% end

for i = 1:22
    switch i
        case  1, color{i} = [233, 40, 52]/255;
        case  2, color{i} = [ 17,168, 97]/255;
        case  3, color{i} = [ 51, 57,146]/255;
        case  4, color{i} = [ 29,185,235]/255;
        case  5, color{i} = [155,165,209]/255;
        case  6, color{i} = [135, 61,147]/255;
        case  7, color{i} = [  0,114,189]/255;
        case  8, color{i} = [237,177, 32]/255;
        case  9, color{i} = [217, 83, 25]/255;
        case 10, color{i} = [126, 47,142]/255;
        case 11, color{i} = [119,172, 48]/255;
        case 12, color{i} = [ 77,190,238]/255;
        case 13, color{i} = [162, 20, 47]/255;
        case 14, color{i} = [  0,114,189]/255;
        case 15, color{i} = [217, 83, 25]/255;
        case 16, color{i} = [237,177, 32]/255;
        case 17, color{i} = [232, 36,142]/255;
        case 18, color{i} = [ 80,187, 79]/255;
        case 19, color{i} = [ 96, 84,168]/255;
        case 20, color{i} = [ 92,198,188]/255;
        case 21, color{i} = [  0,  0,  0]/255;
        case 22, color{i} = [ 52, 62,148]/255;
    end
end

if str2double(erase(version('-release'),{'a', 'b'})) > 2018
    color = cell2mat(color');
    color(color * 1.3 < 1) = color(color * 1.3 < 1) * 1.3;
    colororder(color)
end

if N <= 10
    return;
end

temp = lines(N);

for i = 11:N
    color{i} = temp(i,:);
end

