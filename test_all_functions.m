function test_all_functions()
% 综合测试项目中的多个函数

    % 测试 fitness 函数
    fprintf('开始测试 fitness 函数...\n');
    test_fitness_function();
    
    % 测试 fitness2 函数
    fprintf('\n开始测试 fitness2 函数...\n');
    test_fitness2_function();
    
    % 测试 validate 函数
    fprintf('\n开始测试 validate 函数...\n');
    test_validate_function();
    
    % 测试其他函数（如存在）
    fprintf('\n所有函数测试完成。\n');
end

function test_fitness_function()
% 测试 fitness 函数
    if ~exist('t.mat', 'file') || ~exist('sc.mat', 'file')
        fprintf('  跳过 fitness 函数测试 - 缺少必需的.mat文件\n');
        return;
    end
    
    % 使用合理的参数值进行测试
    test_input = [0.001, 0.0001, 1e-6, 1, 1, -3.1];
    try
        result = fitness(test_input);
        if ~isnan(result) && ~isinf(result)
            fprintf('  fitness 函数测试通过 - 输入: [0.001, 0.0001, 1e-6, 1, 1, -3.1], 输出: %f\n', result);
        else
            fprintf('  fitness 函数测试警告 - 输出为 NaN 或 Inf: %f\n', result);
        end
    catch ME
        fprintf('  fitness 函数测试失败 - 错误: %s\n', ME.message);
    end
end

function test_fitness2_function()
% 测试 fitness2 函数
    if ~exist('t.mat', 'file') || ~exist('sc.mat', 'file')
        fprintf('  跳过 fitness2 函数测试 - 缺少必需的.mat文件\n');
        return;
    end
    
    % 准备 fitness2 函数需要的参数
    Q = [0.020868421 0.016911111 0.010409722 0.016051389 0.020252778];
    l = 60;  % 滤水管上部深度
    d = 72;  % 滤水管下部深度
    M = 74;  % 含水层厚度
    
    test_input = [0.001, 0.0001, 1e-6, 1, 1, -3.1, -3.2, -3.3, -3.4, -3.5];
    try
        result = fitness2(test_input, Q, l, d, M);
        fprintf('  fitness2 函数测试 - 输入: 10维向量, 输出: %s (size: [%d,%d])\n', ...
            class(result), size(result,1), size(result,2));
    catch ME
        fprintf('  fitness2 函数测试失败 - 错误: %s\n', ME.message);
    end
end

function test_validate_function()
% 测试 validate 函数
    if exist('validate.m', 'file')
        if ~exist('t.mat', 'file') || ~exist('sc.mat', 'file')
            fprintf('  跳过 validate 函数测试 - 缺少必需的.mat文件\n');
            return;
        end
        
        % validate函数需要一个参数向量，根据代码应至少包含10个元素
        test_input = [0.0005, 0.0001, 0.0054, 1.4, 0.4, -2.43, -1.58, -1.58, -1.65, -1.71];
        try
            result = validate(test_input);
            fprintf('  validate 函数测试 - 输入: 10维向量, 输出: %f\n', result);
        catch ME
            fprintf('  validate 函数测试失败 - 错误: %s\n', ME.message);
        end
    else
        fprintf('  validate 函数不存在，跳过测试\n');
    end
end