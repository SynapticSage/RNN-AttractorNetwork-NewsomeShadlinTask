function choice = mapChoice(trialvec)
% Simple small function that maps trial input the the desired choice

choice = 'NULL';
if trialvec(1) < 0
    if trialvec(2) > 0
        choice = 'right';
    elseif trialvec(2) < 0
        choice = 'left';
    end
elseif trialvec(1) > 0
    if trialvec(3) > 0
        choice = 'right';
    elseif trialvec(3) < 0
        choice = 'left';
    end
end

end