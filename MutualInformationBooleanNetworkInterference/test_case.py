K = 3

observed_state_sequence = []
observed_state_sequence += [[1, 0, 1, 1, 0, 0, 1, 0, 0, 0]] # t = 1
observed_state_sequence += [[1, 0, 1, 1, 0, 0, 1, 0, 0, 0]] # t = 2
observed_state_sequence += [[1, 0, 1, 1, 0, 0, 1, 0, 0, 0]] # t = 3
observed_state_sequence += [[1, 0, 1, 1, 0, 0, 0, 0, 0, 1]] # t = 4
observed_state_sequence += [[1, 1, 1, 1, 0, 1, 0, 1, 0, 1]] # t = 5
observed_state_sequence += [[0, 1, 0, 1, 1, 1, 0, 1, 1, 1]] # t = 6
observed_state_sequence += [[0, 1, 0, 1, 1, 1, 0, 1, 1, 1]] # t = 7
observed_state_sequence += [[0, 1, 0, 1, 1, 1, 0, 1, 1, 1]] # t = 8
observed_state_sequence += [[0, 1, 0, 1, 1, 1, 0, 1, 1, 1]] # t = 9
observed_state_sequence += [[0, 1, 0, 1, 1, 1, 0, 1, 1, 1]] # t = 10
observed_state_sequence += [[0, 1, 0, 1, 1, 1, 0, 1, 1, 1]] # t = 11
observed_state_sequence += [[0, 1, 0, 0, 1, 1, 0, 1, 1, 1]] # t = 12
observed_state_sequence += [[0, 1, 0, 0, 1, 1, 0, 1, 1, 1]] # t = 13
observed_state_sequence += [[0, 1, 0, 0, 1, 1, 0, 1, 1, 1]] # t = 14
observed_state_sequence += [[0, 1, 0, 0, 1, 1, 0, 1, 1, 1]] # t = 15
observed_state_sequence += [[0, 1, 0, 0, 1, 1, 0, 1, 1, 1]] # t = 16
observed_state_sequence += [[0, 1, 0, 0, 1, 0, 0, 1, 1, 1]] # t = 17
observed_state_sequence += [[0, 1, 0, 0, 1, 0, 0, 1, 1, 1]] # t = 18
observed_state_sequence += [[0, 1, 0, 0, 1, 0, 0, 1, 1, 1]] # t = 19
observed_state_sequence += [[0, 1, 0, 0, 1, 0, 0, 1, 1, 1]] # t = 20
observed_state_sequence += [[0, 1, 0, 0, 1, 0, 0, 1, 1, 1]] # t = 21

observed_state_sequence = [list(x) for x in zip(*observed_state_sequence)]
print(observed_state_sequence)

