# Introduction to shell (bash) scripting

## Sreenu Vattipally
### MRC-University of Glasgow Centre for Virus Research
### University of Glasgow, G61 3RX

---

## What is shell scripting?

Shell scripting refers to writing programs using one of SHELL languages, i,e, Bourne Again SHell (bash), C Shell (ksh), Korn Shell (ksh) or Z shell (zsh). In this course, you will be learning bash scripting. 

Shell scripting offers several benefits, including automation, efficiency, consistency, and flexibility. Here are some of the key advantages:

1. Automation: Bash scripts can automate repetitive tasks, saving time and reducing manual effort. This makes it an ideal tool for system administrators who manage large server infrastructures or users who perform similar actions frequently.

2. Efficiency: Scripts can be run quickly, eliminating the need to manually execute commands multiple times. They also allow you to store and reuse code snippets, further improving efficiency and productivity.

3. Consistency: Bash scripts ensure consistency across systems by allowing for reproducible results regardless of who runs the script or on which machine it is executed. This is particularly useful in environments where multiple people may be performing similar tasks but with different levels of expertise or familiarity with the system.

4. Flexibility: Scripts can be easily modified to accommodate changes in your environment or workflow. 

5. Reduce errors: Scripts will eliminate typographical errors and save debugging time. 



Here's an example of a simple bash script:

```bash
#!/bin/bash
echo "Hello World!"
```

This script prints out the message "Hello World!" when executed. In this, the first line `#!/bin/bash` is the header line (or shebang line), it tells which program to use for executing the following code.  The `echo` command prints the provided text.

Here is an example of printing a variable. 
```bash
#!/bin/bash

Name="Sreenu"
echo "Hello $Name"
```


Here the variable `Name` stores the values `Sreenu`. 
There are several ways to declare variables in BASH scripting:

1. Simple Variable Declaration: This is the most basic way of declaring a variable. You can assign a value directly using the "=" sign like this:

```bash
variable=value
```

2. Command line values to variable
```bash
#!/bin/bash

Name=$1
echo "Hello $Name"
```

3. Command output: You can use command's output values to variables. For example, if you want to set a variable called 'result' to the output of a command, you can do it like this:

```bash
result=$(command)
```

```bash
#!/bin/bash

currentDir=pwd
echo "$currentDir"
currentDir=$(pwd)
echo "$currentDir"
```

(What is the difference between the above  two outputs)

4. Here Document: The "here document" construct is used for setting multi-line strings as the value of a variable. This is done by enclosing the string within single or double quotes and using a heredoc syntax:

```bash
variable=$(cat <<EOF
This is a multiline string
with newlines
and other characters
EOF)

echo "$variable"
```

5. Storing arithmetic values
```bash
var=$((2 + 2))
echo $var
```

---

## String manipulation
### concatenation

Below script will combine two strings (`a` and `b`) and store it in `c`
```bash
#!/bin/bash

a="Hello"
b="World"
c="$a $b"

echo "$c"
```

### Print the length of the string
```bash
#!/bin/bash

a="Hello"
b="World"
c="$a $b"

echo "${#c}" # This will print 11, the length of "Hello World"
```

### Extracting substring
```bash

a="Hello"
b="World"
c="$a $b"

d=${c:3} # Delete the first 3 char
echo $d

d=${c:3:4} # Delete the first 3 char and print 4 after
echo $d

d=${c: -4} # Print the last 5 char
echo $d
```


### String Replacement

To replace a substring within a string:
```bash
str="Hello world"
str=${str/world/universe}
echo $str
```

you can use the `sed` command. For example:
   ```bash
   str="Hello World"
   echo $str | sed 's/World/Universe/' # Output: Hello Universe
   ```

### Changing the case

```bash

string="hello world"

# Convert to uppercase
uppercase=${string^^}
echo $uppercase  # Output: HELLO WORLD


# Convert to lower case
string="Hello World"
lowercase=${string,,}
echo $lowercase;
```

### Remove prefix or suffix
```bash
# To remove preffix use "#"
str=${str#preffix}

# To remove suffix use "%"
str=${str%suffix}

```


## Conditions

Conditions will execute commands in a controlled manner. You can execute commands only if you meet a previous condidtion

```bash
#!/bin/bash

Name=$1

if [ $Name == "Sreenu" ]; then
	echo "Hello $Name"
fi
```

Using else
```bash
#!/bin/bash

Name=$1

if [ $Name == "Sreenu" ]; then
	echo "Hello $Name"
else
	echo "Hello Stranger"
fi
```

Conditions for checking arithmetic values
```bash

var=$1

if [ $var -gt 5 ]; then

	echo "It is greater than 5"

else
	echo "It is less than 5"
fi
```


Here is the list of conditional expressions to check the strings

```bash
$a < $b # Less than
$a > $b # Greater than
$a >= $b # Greater than or equal to
$a <= $b # Less than or equal to
$a = $b # Equal to
$a != $b # Not equal to
```

Here is the list of conditional expressions to check the numerical values

```bash
$a -lt $b # Less than
$a -gt $b # Greater than
$a -ge $b # Greater than or equal to
$a -le $b # Less than or equal to
$a -eq $b # Equal to
$a -ne # Not equal to
num1%num2 # modulo operator – gives the remainder
```

The last expression is called modulo operator. It gives the reminder value of num1 divided with num2


## Loops

Loops are used for running repetitive jobs. There are three different loops in bash scripting.

1. while
2. until
3. for


### Print up to 10 using a while loop
```bash
i=0

while [ $i -le 10 ]
do
	echo i is $i;
	((i++)); # or let i++
done
```

### Read a file line by line and print

```bash

while read line
do
	echo $line
done < $1
```

### Print up to 10 using a until loop

```bash
i=0

while [ $i -ge 10 ]
do
	echo i is $i;
	((i++));
done
```


### Print up to 10 using a for loop

```bash
for ((i=1; i<=10; i++))
do
	echo i is $i;
done
```

### Print up to 10 using a for loop (using range)

```bash
for i in {1..10}
do
	echo i is $i;
done
```


### Print up to 10 using a for loop (using seq)

```bash
for i in $(seq 1 10)
do
	echo i is $i;
done
```


### Print from an array

```bash
array=("element1" "element2" "element3")

for i in "${array[@]}"
do
	echo $i
done
```

### Print directory contents (emulate `ls` command)

```bash
for file in *
do
	echo $file
done
```

### Print specific files (files with “fq” extension)

```bash

for file in *.fq
do
	echo $file
done
```


---

## Exercises
1. Create a reverse complementary genome sequence of `/home/manager/Sreenu/exFiles/Ebola.fa`

2. Print the genome composition (number of different nucleotides) and length of the above genome

3. Convert from fastq to fasta sequence from `/home/manager/Sreenu/FQs/goodData-1.fq`

4. Convert from fastq to fasta sequence only if they have EcoRI (GAATTC) site in the above reads
