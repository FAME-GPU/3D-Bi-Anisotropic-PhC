clc;clear;
fileList = {
    '../test_data_store/woodpile_pcg_n9_gamma_choose4_=0.2000000.mat';
    '../test_data_store/woodpile_pcg_n18_gamma_choose4_=0.2000000.mat';
    '../test_data_store/woodpile_pcg_n36_gamma_choose4_=0.2000000.mat';
    '../test_data_store/woodpile_pcg_n72_gamma_choose4_=0.2000000.mat';
    '../test_data_store/woodpile_pcg_n144_gamma_choose4_=0.2000000.mat';
};

ew_converge = zeros(5, 6);

for i = 1:length(fileList)
    currentFile = fileList{i};

    if ~exist(currentFile, 'file')
        error('文件不存在：%s，请确认文件路径和名称是否正确！', currentFile);
    end

    load(currentFile);
    
    ew_processed = gather(Freq_mtx);
    
    if length(ew_processed) < 6
        error('第%d个文件%s中的ew仅包含%d个特征值，不足6个！', i, currentFile, length(ew_processed));
    end
    
    ew_converge(i, :) = ew_processed(1:6);
end

save('ew_converge_result.mat', 'ew_converge');

fprintf('数据处理完成！\n5行6列矩阵ew_converge已保存为当前目录的ew_converge_result.mat文件\n');
fprintf('矩阵预览：\n');
disp(ew_converge);