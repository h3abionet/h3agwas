
#  Instructions for using AWS Batch

These are supplementary instructions for setting up AWS Batch in order
to use our workflow. These instructions assume that you understand the
basics of AWS.

In summary you must
* do the set-up (this need only be done once -- unless your data set sizes vary very considerably or you change region)
  * decide which AWS region you want to run in
  * create an AWS S3 bucket in that region as work directory (I'll refer to this as `eg-aws-bucket` in these instructions, which    you should replace with your own bucket name
    * if your input and output are going to AWS S3 they should be in the same region -- the input could be on your local disk or S3, and the output could go your local disk or S3. It would probably make sense for the input to be on S3 since you might want to run several times, but that is completely up to you and YMMV. Do remember to set the security on the buckets so as to protect any sensitive data
  * set up an AWS Batch Compute Environment with appropriate permissions
  * set up an AWS Batch queue that connects to the batch environment
  * set up a Nextflow config file that puts all this together
* run the workflow

*NB:* AWS Batch requires all resources to be set up in the same AWS region. A major cause of misconfiguration is failing to observe this (easy to do), so pay attention to which region you are in. In particular, the S3 buckets and AWS Batch needs to be in the same region.



## 1. AWS region and set up

Log into your AWS console.

Choose the AWS region based on what is convenient for you and has a the best pricing. It probably doesn't make a big difference. Our University's ethics committee prefers us to run on `af-south-1` since these are subject to South African privacy laws. Remember that any buckets you use must be in the same region.

You need to have a bucket that can be used for temporary storage. It could be an existing bucket. Make sure you are the only person who can read and write to it.






## 2. Set up your Batch Compute Environment

This defines what resources you need for your workflow.  There are a number of steps but mainly we can use the default values. However, you do need to be careful that you set things up so that you have enough reources.

### 2.1 Disk space

By default, AWS instances that run batch jobs are 30GB in size. We think that you need an image that is at least 4x bigger than the input data size to run   safely. If your need less than 30GB, there's no problem _and you can skip the rest of 2.1_. If not, there's an extra configuration step in the configuration to set up an environment with disks of the correct size

The easiest way of doing this is to set up a _launch template_. You can do this using the console but in my experience it is more complex than using  the command line tools.

Install boto3 library using pip or yum or the like (e.g., `yum install python3-boto`)

There is a file called `launch-template.json`  in this directory. Download it to your machine. Change the _LaunchTemplateName_ field associated value to something meaningful and unique for you and the _VolumeSize_ field to the value you want (in units of GB)

Then, using the following command (_mutatis mutandis_ -- that is change the region and the name of the template file if you've changed it) create the launch template


```
 aws ec2 --region af-south-1 create-launch-template --cli-input-json file://h3a-100GB-template.json
```


### 2.2. RAM size

The bigger your data files the more RAM you need. This is not an issue you have to worry about at configuration time. When you *run* the workflow you may have to set the `plink_mem_req` and `other_mem_req` parameters as discussed elsewhere.

## 2.3 Setting up the environment and queue

Read

*  https://www.nextflow.io/docs/latest/executor.html#aws-batch
*  https://www.nextflow.io/docs/latest/awscloud.html#awscloud-batch  (up to point 5 -- from "Configuration" on is meant for developers of pipelines and not users).


### 2.3.1 Create a role that set up permissions

Go to IAM to set up a policy that allows the AWS service to call Batch on your behalf and to access Amazon S3
* Pick Roles
* Choose Create Roles
* Pick S3 and then the S3 Batch Operations use case
* Next: Permissions
  * Type in AmazonS3FullAccess in the filter and tick
  * Choose Next:Tags
* Choose Next: Review
  * Give a meaninful name (e.g. GWASBatchRole)
  * Choose "Create Role"
* Then click on that role to see the Summary
  * Choose _Attach policies_
  * In the filter, choose "AWSBatchServiceRole" and select it (tick it)
  * Choose _Attach Policy"


### 2.3.1 Setting up an environment
You should be able to follow the default settings unless you need a launch template

Choose
* Computer Environment Configuration
   * _Managed_ environment type
   * give a meaningful name
   * enable environment
   * under additional settings
        * _AWSBatchServiceRole_
	* ecsInstanceRole
	* Choose a keypair for the region (not usually needed)
* Instance Configuration
  * _Spot_ pricing (choose the percentage your pocket can afford)
  * _Optimal_ for allowed instance types
  * `SPOT_CAPACITY_OPTIMIZED` for allocation strategy
  * Under _Additional settings_
     * pick _none_ or a template you have defined if you need to. Note for each template you need to define a new environment. (See section 2.1 above)
    * You don't need to to pick an AMI and should do so only if you really know what you are doing.
* Under networking you can pick the defaults -- note that under additonal settings are the definitions of the securty groups which define access
* Add tags if you want to  : may be helpful for tracking resources


### 2.3.2 Adding S3 access permissions

The AWS Batch instances will need to access S3 and so you need to give them this permission. 
* From the AWS Console, choose "Services" and then "IAM"
* Choose "Roles"
* Choose "ecsInstanceRole"
* Choose `Attach Policies`
* In the filter bar type in `AmazonS3FullAccess` (NB: no spaces) and select

The _ecInstanceRole_ now comprises two policies: _AmazonEC2ContainerServiceforEC2Role_ and _AmazonS3FullAccess_


### 2.3.3 Set up a queue

In the Amazon Console, go back to _AWS Batch_ from the list of services

Create a job queue. Unless you need to do something fancy just pick the default options.
* for convenience call the queue the same as the environment
* attach the environment you created to the queue

You use the queue name in your Nextflow config file.


## 3. Create an additonal Nextflow configuration file

It will look something like this. I'll call this `aws.config` for the example but you can name it what ever you want.

```

params {
   plink_mem_req = "8GB"
   other_mem_req = "8GB"
}


profiles {


    awsbatch {
         region = "af-south-1"
         accessKey ='YourAccessKeyForTheRegion'
         secretKey ='AssociatedSecretKey'
         aws.uploadStorageClass = 'ONEZONE_IA'
         process.queue = 'QueueYouCreated'
         process.executor = "awsbatch"

    }


}
```


Then run the code something like this. Note in this example I am using input  that's already in S3, except for the one file is local (this is to show you that the data can be in different places -- it would probably make
sense for the phenotype file to stored in S3, but perhaps you are trying to be extra careful).

```
nextflow run h3agwas/qc/main.nf \
       -profile awsbatch -c aws.config \
       -work-dir 's3://za-batch-work' \
       --input_dir 's3://za-batch-data' --input_pat sim1_s \
       --output_dir 's3://za-batch-data/qc' --output sim1_s_qc\
       --data 's3://za-batch-data/out_qt.pheno' \
       --case_control data/see.fam --case_control_col=pheno
```

## 4 Clean up

Note the work bucket you give will start to fill up (as will any output buckets). If you do lots of analysis it's possible for the work bucket to quickly get to the hundreds of GB level. There may be some sensitive data and AWS will also charge you for this, so remember to regularly delete objects from your work bucket.
