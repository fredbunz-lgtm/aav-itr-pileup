import os
import requests
from dotenv import load_dotenv

load_dotenv()

token = os.environ.get("ADDGENE_API_TOKEN")

if not token:
    raise RuntimeError("ADDGENE_API_TOKEN not found in .env")

base_url = "https://www.addgene.org/api/"
headers = {
    "Authorization": f"Token {token}",
    "Accept": "application/json",
}

response = requests.get(base_url, headers=headers, timeout=30)

print("Status code:", response.status_code)
print("Content-Type:", response.headers.get("Content-Type"))

try:
    data = response.json()
    print("JSON response keys:", list(data) if isinstance(data, dict) else type(data))
except Exception:
    print(response.text[:1000])
