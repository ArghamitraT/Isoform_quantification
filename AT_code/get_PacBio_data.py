from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
import time
import requests
import os


def download_file(url, local_filename):
    """Download a file from a URL into the specified directory."""
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    #return local_filename

def load_url(driver, url):

    # Navigate to the page
    driver.get(url)

    # Wait for the page to load dynamically; adjust time as necessary
    time.sleep(10)  # Sleep to allow page contents to load

    # Find all <a> elements and filter for .bam files using By strategy
    elements = driver.find_elements(By.TAG_NAME, 'a')

    file_urls = [element.get_attribute('href') for element in elements if
                 element.get_attribute('href')]

    return file_urls
def get_file_urls(url):
    # Setup Selenium WebDriver
    service = Service(ChromeDriverManager().install())
    options = webdriver.ChromeOptions()
    driver = webdriver.Chrome(service=service, options=options)

    # Getting the folders for each day
    file_urls =  load_url(driver, url)

    # filtering out the right folders
    right_url = []
    for url in file_urls:
        if url.split('/')[-1] == '':
             if url.split('/')[-3].split('-')[0][:3] == 'day' and url.split('/')[-3].split('-')[1][:3] == 'rep':
                 right_url.append(url)
                 local_dir = home_dir + '/'.join(url.split('/')[-3:-1])
                 if not os.path.exists(local_dir):
                     os.makedirs(local_dir)


    # getting the url inside day-rep folder
    for bam_url in right_url:
        file_urls = load_url(driver, bam_url)

        # filtering out the right folders
        for url_F in file_urls[5:]:
            local_filename = home_dir + '/'.join(url_F.split('/')[-3:])
            download_file(url_F, local_filename)

        print()


    driver.quit()
    return file_urls


PacBio_url = 'https://downloads.pacbcloud.com/public/dataset/MAS-Seq-bulk-2023-WTC11/'
home_dir = '/Users/arghamitratalukder/Google Drive/My Drive/David_Argha_shared/RNA_Splicing/data/PacBio_data_Liz/'


PacBio_days = 'day0-rep1/'
# The public URL
bam_file_urls = get_file_urls(url=PacBio_url+PacBio_days)
for url in bam_file_urls:
    print(f"Downloading {url}...")
    download_file(url, home_directory)
    print(f"Saved to {home_directory}")

print()

# bam_file_urls = [element.get_attribute('href') for element in elements if
#                  element.get_attribute('href') and element.get_attribute('href').endswith('.bam')]

