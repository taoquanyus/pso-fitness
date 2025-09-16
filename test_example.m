function test_example()
    % 这是一个示例测试用例文件
    % 用于测试 basic addition 功能

    % 测试用例 1: 验证 1 + 1 = 2
    result = 1 + 1;
    expected = 2;
    assert(result == expected, 'Test 1 Failed: 1 + 1 should equal 2');

    % 测试用例 2: 验证 2 + 3 = 5
    result = 2 + 3;
    expected = 5;
    assert(result == expected, 'Test 2 Failed: 2 + 3 should equal 5');

    fprintf('All tests passed!\n');
end