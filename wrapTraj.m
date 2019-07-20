function wrappedTraj = wrapTraj(traj)

% Find jumps
posJumpIdx = find(diff(traj)>100);

for i = 1:numel(posJumpIdx)
    if i == 1
        traj(1:posJumpIdx(i)) = traj(1:posJumpIdx(i)) + 180;
    else
        traj(posJumpIdx(i-1):posJumpIdx(i)) = traj(posJumpIdx(i-1):posJumpIdx(i)) + 180;
    end
end


negJumpIdx = find(diff(traj)<-100);

for i = 1:numel(negJumpIdx)
    if i == 1
        traj(1:negJumpIdx(i)) = traj(1:negJumpIdx(i)-1) - 180;
    else
        traj(negJumpIdx(i-1):negJumpIdx(i)) = traj(negJumpIdx(i-1):negJumpIdx(i)) - 180;
    end
end


wrappedTraj = traj;
end