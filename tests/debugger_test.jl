x = 1
y = 2
z = x+y
println("z = "*string(z))

function ftest(x)
    println("Entering ftest")
    for i in 1:10
        x += 1
    end
    print("Done")
end

print("Test of debugging functions")
@enter y = ftest(5)