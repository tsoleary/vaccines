 function[state,options,optchanged] = outputFcn(options, state, flag)
% %  optchanged = false;
% %  disp(state.LastImprovement)
% % 
% %  Output.LastImprovement=state.LastImprovement;
  persistent state_record 
    if isempty(state_record)
      state_record = struct('Population', {}, 'Best', {}, 'Score', {}, 'LastImprovement', {});
    end
    if nargin == 0
      state = state_record;
      options = [];
      optchanged = [];
    else
      state_record(end+1) = struct('Population', state.Population, 'Best', state.Best', 'Score', state.Score, 'LastImprovement', state.LastImprovement);
      optchanged = false;
    end