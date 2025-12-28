% MergeRefactoredCode.m
% 功能：将当前目录下的 stapmat.m 以及 SRC 文件夹及其子文件夹中的所有 .m 文件合并
% 目的：用于全代码审查
% 输出文件：All_MATLAB_Full_Project.txt

clear; clc;

% --- 配置 ---
targetFolder = 'SRC'; % 类库文件夹
mainFileName = 'stapmat.m'; % 主函数文件名
outputFileName = 'All_MATLAB_Full_Project.txt'; % 输出文件名

% --- 开始 ---
fprintf('正在准备整合项目代码...\n');

% 1. 获取主函数信息
mainFile = dir(mainFileName);

% 2. 递归获取 Classes 文件夹下的所有 .m 文件
classFiles = dir(fullfile(targetFolder, '**', '*.m'));

% 3. 合并文件列表 (主函数排在最前面)
files = [mainFile; classFiles];

if isempty(files)
    fprintf('错误: 未找到任何代码文件。\n');
    return;
end

% 打开输出文件
fidOut = fopen(outputFileName, 'w');
if fidOut == -1
    error('无法创建输出文件: %s', outputFileName);
end

fprintf('找到 %d 个文件，开始合并...\n', length(files));

for i = 1:length(files)
    % 获取完整路径
    fullPath = fullfile(files(i).folder, files(i).name);
    
    % 计算显示路径
    % 如果是类文件，显示相对于项目根目录的路径；如果是主函数，仅显示文件名
    k = strfind(fullPath, targetFolder);
    if ~isempty(k)
        relativePath = fullPath(k(1):end);
    else
        relativePath = files(i).name;
    end
    
    % 读取内容
    fidIn = fopen(fullPath, 'r');
    if fidIn == -1
        fprintf('警告: 无法读取文件 %s，跳过。\n', fullPath);
        continue;
    end
    fileContent = fread(fidIn, '*char')';
    fclose(fidIn);
    
    % 写入分隔符和文件名头
    fprintf(fidOut, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(fidOut, '%% FILE PATH: %s\n', relativePath);
    fprintf(fidOut, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    
    % 写入内容
    fprintf(fidOut, '%s\n\n', fileContent);
    
    fprintf('已合并: %s\n', relativePath);
end

% 关闭输出文件
fclose(fidOut);

fprintf('\n合并完成！\n');
fprintf('输出文件已生成: %s\n', fullfile(pwd, outputFileName));
fprintf('您可以将此文件发送给我进行一致性审查。\n');