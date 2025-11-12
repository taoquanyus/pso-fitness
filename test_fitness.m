function test_fitness()
% 测试 fitness 函数的测试用例
% 注意：fitness函数依赖于t.mat和sc.mat文件，这些文件必须存在

    % 检查必需的.mat文件是否存在
    if ~exist('t.mat', 'file') || ~exist('sc.mat', 'file')
        fprintf('错误：缺少必需的t.mat或sc.mat文件\n');
        return;
    end
    
    % 根据fitness函数的实现，输入x应包含至少5个元素 (x(1)到x(5))
    % x(1) = T/M, x(2) = S (储水系数), x(3) = a (导水系数), x(4) = P, x(5) = C
    % 还需要x(6)作为alpha(3)的值
    
    % 测试用例1：基本功能测试
    test_input1 = [0.001, 0.0001, 1e-6, 1, 1, -3.1];
    try
        result1 = fitness(test_input1);
        fprintf('测试用例1 - 输入: [0.001, 0.0001, 1e-6, 1, 1, -3.1], 输出: %f\n', result1);
    catch ME
        fprintf('测试用例1 - 错误: %s\n', ME.message);
    end
    
    % 测试用例2：不同参数值测试
    test_input2 = [0.01, 0.001, 1e-5, 2, 0.5, -4.0];
    try
        result2 = fitness(test_input2);
        fprintf('测试用例2 - 输入: [0.01, 0.001, 1e-5, 2, 0.5, -4.0], 输出: %f\n', result2);
    catch ME
        fprintf('测试用例2 - 错误: %s\n', ME.message);
    end
    
    % 测试用例3：小数值测试
    test_input3 = [1e-5, 1e-6, 1e-8, 0.5, 0.1, -2.5];
    try
        result3 = fitness(test_input3);
        fprintf('测试用例3 - 输入: [1e-5, 1e-6, 1e-8, 0.5, 0.1, -2.5], 输出: %f\n', result3);
    catch ME
        fprintf('测试用例3 - 错误: %s\n', ME.message);
    end
    
    % 测试用例4：边界值测试
    test_input4 = [1, 1, 1, 1, 1, -3.2];
    try
        result4 = fitness(test_input4);
        fprintf('测试用例4 - 输入: [1, 1, 1, 1, 1, -3.2], 输出: %f\n', result4);
    catch ME
        fprintf('测试用例4 - 错误: %s\n', ME.message);
    end
    
    fprintf('所有测试用例执行完毕。\n');
end